import FlowCytometryTools as FCT
import json,os
from synthetic_circuit_performance import *
from rank_order_truth_tables import *
import numpy as np
import itertools
import pandas as pd
import ast,random
from pprint import pprint

input_states=["00","01","10","11"]

def get_channels(lab_name):
    if "BIOFAB" in lab_name:
        GFP = "FL1-A"
        Sytox = "FL4-A"
        forward_scatter = "FSC-A"
    elif "Transcriptic" in lab_name:
        GFP = "BL1"
        Sytox = "RL1"
        forward_scatter = "FSC"
    else:
        raise ValueError("Lab name {} not implemented.".format(lab_name))
    return {"GFP" : GFP,"Sytox" : Sytox,"FSC" : forward_scatter }


def format_data(prefix,fname,transform=None,channels=None):
    sample = FCT.FCMeasurement(ID="temp",datafile=os.path.expanduser(os.path.join(prefix,fname)))
    if transform == 'hlog':
        # see FlowCytometryTools documentation
        tsample = sample.transform(transform,channels=[c for _,c in channels.items()], b=100)
    elif transform == 'tlog':
        tsample = sample.transform(transform,channels=[c for _,c in channels.items()], th=2)
    elif transform is not None:
        raise ValueError("Transform {} not recognized.".format(transform))
    return tsample if transform else sample
    

def get_data_tx(circuit,ingest_file="transcriptic_april_fcsfiles_dan.csv",transform=None):
    df = pd.read_csv(open(ingest_file))
    data = {}
    for index, row in df.iterrows():
        if row["bead"]:
            continue
        elif getcircuit(row["gate"]) == circuit:
            lab_name = "Transcriptic"
            channels = get_channels(lab_name)
            fname = row["fcs_files"][0]
            media = row["media"].split('/')[-2]
            input_state = "".join([str(s) for s in ast.literal_eval(row["input"])])
            rep = row["replicate"]
            od = row["od"]
            metadata = {'media':media,'circuit':circuit,'od':od,'input_state':input_state,'rep':rep}
            sample = format_data("",fname,transform,channels)
            d = fname.split('/')[-4]
            if d in data:
                data[d].append((sample,channels,metadata)) 
            else:
                data[d] = [(sample,channels,metadata)]  
    return data


def get_data_mongo(ingest_file, prefix,transform=None):
    matches = json.load(open(ingest_file))
    data = {}
    for m in matches:  
        lab_name = m['attributes']['lab']
        channels = get_channels(lab_name)
        fname = m['measurement']['file']['filename']
        media = m['contents'][0]['name']['label']
        circuit = m['strain']['circuit']
        input_state = m['strain']['input_state']
        try:
            rep = m['replicate']
        except:
            rep = None
        try:
            od = m['inoculation_density']['value']
        except:
            od = None
        metadata = {'media':media,'circuit':circuit,'od':od,'input_state':input_state,'rep':rep}
        sample = format_data(prefix,fname,transform,channels)
        d = fname.split('/')[-2]
        if d in data:
            data[d].append((sample,channels,metadata)) 
        else:
            data[d] = [(sample,channels,metadata)]  
    return data
            

def get_data(mongo=True,ingest_file= "matches_biofab_nand.json", prefix="~/sd2e-community/uploads",transform=None,circuit=None):
    # transform = 'hlog'
    if mongo:
        data = get_data_mongo(ingest_file, prefix, transform)
    else:
        data = get_data_tx(circuit,transform=transform)
    return data


def bin_data(pts,bin_endpoints):
    hist = [0]*(len(bin_endpoints)+1)
    for p in list(pts):
        ind = sorted(bin_endpoints+[p]).index(p)
        hist[ind]+=1
    return hist

def get_bin_centers(bin_endpoints):
    return np.asarray([a + (b-a)/2 for (a,b) in zip([0]+bin_endpoints,bin_endpoints+[2*bin_endpoints[-1] - bin_endpoints[-2]])])
        

def sort_strains_into_histograms(data,bin_endpoints=list(range(250,10250,250)),media="",od=""):
    circuits = {}
    for key,samp_list in data.items():
        for sample,channels,metadata in samp_list:
            if (not media or metadata["media"] == media) and (not od or metadata['od']==od):
                print(metadata)
                ip = metadata["input_state"]
                md = metadata.copy()
                md.pop("input_state")
                md = tuple(md.items())
                pts = list(sample.data[channels["GFP"]].values.transpose())
                hist = bin_data(pts,bin_endpoints)
                if md not in circuits:
                    circuits[md]  = {"00" : [], "01" : [], "10" : [], "11" : []}
                circuits[md][ip].append(hist)
    bin_centers = get_bin_centers(bin_endpoints)
    return circuits, bin_centers


def calculate_results_and_null(data,bin_endpoints=[np.log10(r) for r in range(250,10250,250)],media="Synthetic_Complete",od=0.0003,dump=True):
    circuits, bin_centers = sort_strains_into_histograms(data,bin_endpoints=bin_endpoints,media=media,od=od)
    if dump:
        circuitstr = {str(k):circuits[k] for k in circuits}
        json.dump(circuitstr,open("test_circuits.json","w"))
    results = get_results(circuits,bin_centers)
    null_model = build_null_model(circuits,bin_centers,num_choices=2)
    return results, null_model


def load_circuits(fname):
    circuitstr = json.load(open("test_circuits.json"))
    circuits = {}
    for metadata,vals in circuitstr.items():
        metadata = ast.literal_eval(metadata)
        vals = {k: [np.asarray(hist) for hist in vals[k]] for k in vals}
        circuits[metadata] = vals
    return circuits


def get_results(circuits,bin_centers,num_choices=2):
    new_circuits = {}
    for metadata,vals in circuits.items():
        md = []
        for m in metadata:
            if m[0] not in ["rep","od"]:
                md.append(m)
        md = tuple(md)
        if md not in new_circuits:
            new_circuits[md] = vals
        else:
            temp = dict(new_circuits[md])
            new_circuits[md] = {k: temp.get(k) + vals.get(k) for k in vals.keys()}
    all_scores = { md : {'separations':[],'truthtable_correct':[]} for md in new_circuits}
    for md,ip in new_circuits.items():
        for m in md:
            if m[0] == "circuit":
                desiredtt = desired_truth_tables[m[1]]
                break
        separations = []
        truthtable_correct = []
        # put itertools.product back
        for _ in range(num_choices):
            choice = [random.choice(ip[input_states[k]]) for k in range(4)]
            d = dict(zip(input_states,choice))
            scores = rank_noncst_tables(d,bin_centers)
            separations.append(scores[1][0] - scores[0][0])
            if scores[0][1] == desiredtt:
                truthtable_correct.append(scores[1][0] - scores[0][0])
        all_scores[md]['separations'].extend(separations)
        all_scores[md]['truthtable_correct'].extend(truthtable_correct)
    return all_scores


def build_null_model(circuits,bin_centers,num_choices=2):
    pass


if __name__ == "__main__":
    fname = "test_circuits.json"
    circuits = load_circuits(fname)
    bin_endpoints =[np.log10(r) for r in range(250,10250,250)]
    bin_centers = get_bin_centers(bin_endpoints)
    results = get_results(circuits,bin_centers)
    print(results)












