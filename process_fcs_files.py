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
        GFP = "BL1-A"
        Sytox = "RL1"
        forward_scatter = "FSC-A"
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
        if pd.isnull(row["bead"]) and getcircuit(row["gate"]) == circuit:
            lab_name = "Transcriptic"
            channels = get_channels(lab_name)
            fname = ast.literal_eval(row["fcs_files"])[0]
            fname = os.path.join(os.path.expanduser("~/sd2e-community"),"/".join(fname.split('/')[3:]))
            media = row["media"].split('/')[-2]
            input_state = "".join([str(s) for s in ast.literal_eval(row["input"])])
            rep = row["replicate"]
            od = row["od"]
            metadata = {'media':media,'circuit':circuit,'od':od,'input_state':input_state,'rep':rep}
            channels.pop("Sytox")
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


def get_log_values(values):
    log_vals = []
    for v in values:
        if v > 0:
            log_vals.append(np.log10(v))
        else:
            log_vals.append(-1)
    return np.asarray(log_vals)
            

def get_data(mongo=True,ingest_file= "matches_biofab_nand.json", prefix="~/sd2e-community/uploads",transform=None,circuit=None):
    # transform = 'hlog'
    if mongo:
        data = get_data_mongo(ingest_file, prefix, transform)
    else:
        data = get_data_tx(circuit,ingest_file=ingest_file,transform=transform)
    return data


def bin_data(pts,bin_endpoints):
    hist = [0]*(len(bin_endpoints)+1)
    for p in list(pts):
        ind = sorted(bin_endpoints+[p]).index(p)
        hist[ind]+=1
    return hist

def get_bin_centers(bin_endpoints):
    return np.asarray([a + (b-a)/2 for (a,b) in zip([0]+bin_endpoints,bin_endpoints+[2*bin_endpoints[-1] - bin_endpoints[-2]])])
        

def sort_strains_into_histograms(data,bin_endpoints=list(np.arange(0.05,7.05,0.05)),media=None,od=None):
    circuits = {}
    for key,samp_list in data.items():
        for sample,channels,metadata in samp_list:
            if (not media or metadata["media"] == media) and (not od or metadata['od']==od):
                print(metadata)
                ip = metadata["input_state"]
                md = metadata.copy()
                md.pop("input_state")
                md = tuple(md.items())
                pts = get_log_values(sample.data[channels["GFP"]].values.transpose())
                hist = bin_data(pts,bin_endpoints)
                if md not in circuits:
                    circuits[md]  = {"00" : [], "01" : [], "10" : [], "11" : []}
                circuits[md][ip].append(hist)
    bin_centers = get_bin_centers(bin_endpoints)
    return circuits, bin_centers


def calculate_results_and_null(data,bin_endpoints=[np.log10(r) for r in range(250,10250,250)],media=None,od=None,dumpfile=None,num_choices_res=10,num_choices_null=20):
    print("Initial sort....")
    circuits, bin_centers = sort_strains_into_histograms(data,bin_endpoints=bin_endpoints,media=media,od=od)
    if dumpfile:
        circuitstr = {str(k):circuits[k] for k in circuits}
        json.dump(circuitstr,open(dumpfile,"w"))
    print("Initial sort done.")
    print("Processing results....")
    results = get_results(circuits,bin_centers,num_choices=num_choices_res)
    print("Processing results done.")
    print("Building null model....")
    null_model = build_null_model(circuits,bin_centers,num_choices=num_choices_null)
    print("Building null model done.")
    return results, null_model


def load_circuits(fname):
    circuitstr = json.load(open(fname))
    circuits = {}
    for metadata,vals in circuitstr.items():
        metadata = ast.literal_eval(metadata)
        vals = {k: [np.asarray(hist) for hist in vals[k]] for k in vals}
        circuits[metadata] = vals
    return circuits


def get_results(circuits,bin_centers,num_choices=25):
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
    all_scores = { md : {'truthtable_incorrect':[],'truthtable_correct':[]} for md in new_circuits}
    for md,ip in new_circuits.items():
        for m in md:
            if m[0] == "circuit":
                desiredtt = desired_truth_tables[m[1]]
                break
        truthtable_incorrect = []
        truthtable_correct = []
        for _ in range(num_choices):
            choice = [np.asarray(random.choice(ip[input_states[k]])) for k in range(4)]
            d = dict(zip(input_states,choice))
            scores = rank_noncst_tables(d,bin_centers)
            if scores[0][1] == desiredtt:
                truthtable_correct.append(scores[1][0] - scores[0][0])
            else:
                truthtable_incorrect.append(scores[1][0] - scores[0][0])
        all_scores[md]['truthtable_incorrect'].extend(truthtable_incorrect)
        all_scores[md]['truthtable_correct'].extend(truthtable_correct)
    return all_scores


def build_null_model(circuits,bin_centers,num_choices=50):
    random_circuits = {}
    for metadata,vals in circuits.items():
        md = []
        for m in metadata:
            if m[0] not in ["rep","od"]:
                md.append(m)
        md = tuple(md)
        all_inputs = [ a  for k in vals for a in vals[k] ]
        if md not in random_circuits:
            random_circuits[md] = all_inputs
        else:
            random_circuits[md].extend(all_inputs)
    all_scores = { md : {'truthtable_incorrect':[],'truthtable_correct':[]} for md in random_circuits}
    for md,ip in random_circuits.items():
        for m in md:
            if m[0] == "circuit":
                desiredtt = desired_truth_tables[m[1]]
                break
        truthtable_incorrect = []
        truthtable_correct = []
        for _ in range(num_choices):
            choice = [np.asarray(random.choice(ip)) for _ in range(4)]
            d = dict(zip(input_states,choice))
            scores = rank_noncst_tables(d,bin_centers)
            if scores[0][1] == desiredtt:
                truthtable_correct.append(scores[1][0] - scores[0][0])
            else:
                truthtable_incorrect.append(scores[1][0] - scores[0][0])
        all_scores[md]['truthtable_incorrect'].extend(truthtable_incorrect)
        all_scores[md]['truthtable_correct'].extend(truthtable_correct)
    return all_scores


def test(out=False):
    from plot_fcs_files import plot_bar
    fname = "test_circuits_biofab_SC_od3.json"
    circuits = load_circuits(fname)
    bin_endpoints =[np.log10(r) for r in range(250,10250,250)]
    bin_centers = get_bin_centers(bin_endpoints)
    for md,vals in circuits.items():
        ptclouds = [vals[k][0] for k in input_states]
        colors = ["blue","blue","blue","red"]
        plot_bar(bin_centers,ptclouds,xlim=[1,7],ylim=[0,1],colors=colors)



    results={}
    print("Processing results...")
    results = get_results(circuits,bin_centers)
    ptcloud = results[list(results.keys())[0]]["separations"]
    N = len(ptcloud)
    print(N)
    print("Results calculated.")
    bin_endpts = list(np.arange(0,0.0105,0.0005))
    print("Plotting results...")
    plot_hist_from_point_cloud([ptcloud],xlim=[0,0.01],ylim=[0,len(ptcloud)],bin_endpts=bin_endpts)
    print("Building null model...")
    N = 1296
    nullmod = build_null_model(circuits,bin_centers,num_choices=10*N)
    print("Null model built.")
    print("Plotting null model...")
    ptcloud = nullmod[list(nullmod.keys())[0]]["separations"]
    plot_hist_from_point_cloud([ptcloud], xlim=[0, 0.01], ylim=[0, len(ptcloud)], bin_endpts=bin_endpts)
    if out:
        return results, nullmod

if __name__ == "__main__":
    test()











