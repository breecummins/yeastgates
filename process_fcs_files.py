import FlowCytometryTools as FCT
import json, os
from synthetic_circuit_performance import *
from rank_order_truth_tables import *
import numpy as np
import itertools
import pandas as pd
import ast, random
from pprint import pprint

input_states = ["00", "01", "10", "11"]


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
    return {"GFP": GFP, "Sytox": Sytox, "FSC": forward_scatter}


def format_data(prefix, fname, transform=None, channels=None):
    sample = FCT.FCMeasurement(ID="temp", datafile=os.path.expanduser(os.path.join(prefix, fname)))
    if transform == 'hlog':
        # see FlowCytometryTools documentation
        tsample = sample.transform(transform, channels=[c for _, c in channels.items()], b=100)
    elif transform == 'tlog':
        tsample = sample.transform(transform, channels=[c for _, c in channels.items()], th=2)
    elif transform is not None:
        raise ValueError("Transform {} not recognized.".format(transform))
    return tsample if transform else sample


def make_threshold_gate(sample, threshold, channel, region="above"):
    gate = FCT.ThresholdGate(threshold, [channel], region=region)
    return sample.gate(gate)


def get_data_tx(circuit,ingest_file="transcriptic_april_fcsfiles_dan.csv", prefix="~/sd2e-community", transform=None):
    df = pd.read_csv(open(ingest_file))
    data = {}
    count = 0
    for index, row in df.iterrows():
        if pd.isnull(row["bead"]):
            c = getcircuit(row["gate"])
            if c == circuit:
                lab_name = "Transcriptic"
                channels = get_channels(lab_name)
                fname = ast.literal_eval(row["fcs_files"])[0]
                fname = os.path.join(os.path.expanduser(prefix), "/".join(fname.split('/')[3:]))
                media = row["media"].split('/')[-2]
                try:
                    input_state = "".join([str(s) for s in ast.literal_eval(row["input"])])
                except:
                    # wild type have no input state
                    continue
                rep = row["replicate"]
                od = row["od"]
                metadata = {'media': media, 'circuit': circuit, 'od': od, 'input_state': input_state, 'rep': rep}
                channels.pop("Sytox")
                sample = format_data("", fname, transform, channels)
                d = fname.split('/')[-4]
                if d in data:
                    data[d].append((sample, channels, metadata))
                else:
                    data[d] = [(sample, channels, metadata)]
                count += 1
                if not count % 100:
                    print(count)
    print("Total files chosen = {}".format(count))
    return data


def get_data_mongo(circuit,ingest_file="matches_biofab_all_circuits_all_media.json", prefix="~/sd2e-community/uploads",
                   transform=None):
    matches = json.load(open(ingest_file))
    data = {}
    count = 0
    for m in matches:
        c = getcircuit(m['strain']['circuit'])
        if c == circuit:
            lab_name = m['attributes']['lab']
            channels = get_channels(lab_name)
            fname = m['measurement']['file']['filename']
            media = m['contents'][0]['name']['label']
            input_state = m['strain']['input_state']
            try:
                rep = m['replicate']
            except:
                rep = None
            try:
                od = m['inoculation_density']['value']
            except:
                od = None
            metadata = {'media': media, 'circuit': circuit, 'od': od, 'input_state': input_state, 'rep': rep}
            sample = format_data(prefix, fname, transform, channels)
            d = fname.split('/')[-2]
            if d in data:
                data[d].append((sample, channels, metadata))
            else:
                data[d] = [(sample, channels, metadata)]
            count += 1
            if not count % 100:
                print(count)
    print("Total files chosen = {}".format(count))
    return data


def get_log_values(values):
    log_vals = []
    for v in values:
        if v > 0:
            log_vals.append(np.log10(v))
        else:
            log_vals.append(-1)
    return np.asarray(log_vals)


def get_data(circuit,mongo=True, ingest_file="matches_biofab_nand.json", prefix="~/sd2e-community/uploads",transform=None,threshold=None,channel = "FSC",region="above"):
    # transform = 'hlog'
    def make_thresh(data,threshold):
        if threshold:
            gated_data = {}
            for d,l in data.items():
                new_l = []
                for (sample,channels,metadata) in l:
                    gated_sample = make_threshold_gate(sample,threshold,channels[channel],region=region)
                    new_l.append((gated_sample,channels,metadata))
                gated_data[d] = new_l
            data = gated_data
        return data
    if mongo:
        data = get_data_mongo(circuit,ingest_file, prefix, transform)
    else:
        data = get_data_tx(circuit,ingest_file, prefix, transform)
    if threshold:
        gdata = make_thresh(data,threshold)
        data = gdata
    return data


def bin_data(pts, bin_endpoints):
    hist = [0] * (len(bin_endpoints) + 1)
    for p in list(pts):
        ind = sorted(list(bin_endpoints) + [p]).index(p)                
        hist[ind] += 1
    return hist


def get_bin_centers(bin_endpoints):
    return [a + (b - a) / 2 for (a, b) in zip([0] + list(bin_endpoints), list(bin_endpoints) + [2 * bin_endpoints[-1] - bin_endpoints[-2]])]


def sort_strains_into_histograms(data, bin_endpoints=list(np.arange(0.05, 7.05, 0.05)), media=None, od=None):
    circuits = {}
    for key, samp_list in data.items():
        for sample, channels, metadata in samp_list:
            if (not media or metadata["media"] == media) and (not od or metadata['od'] == od):
                print(metadata)
                ip = metadata["input_state"]
                md = metadata.copy()
                md.pop("input_state")
                md = str(tuple(md.items()))
                pts = get_log_values(sample.data[channels["GFP"]].values.transpose())
                hist = bin_data(pts, bin_endpoints)
                if md not in circuits:
                    circuits[md] = {"00": [], "01": [], "10": [], "11": []}
                circuits[md][ip].append(hist)
    return circuits


def calculate_results_and_null(mongo=True,ingestfile="matches_biofab_all_circuits_all_media.json", prefix="~/sd2e-community/uploads", transform="hlog",bin_endpoints=[np.log10(r) for r in range(250, 10250, 250)], media=None, od=None,dumpfile="temp.json", num_choices_res=10, num_choices_null=20, loadfile=None, circuits=["XNOR","XOR","NOR","OR","NAND","AND"],threshold=None,channel = "FSC",region="above",savefile=None):
    if loadfile:
        cts = json.load(open(loadfile))
        bin_centers = cts["bin_centers"]
        cts.pop("bin_centers")
    else:
        cts = {}
        bin_centers = get_bin_centers(bin_endpoints)
        for circuit in circuits:
            print("Getting data for circuit {}....".format(circuit))
            data = get_data(circuit, mongo=mongo, ingest_file=ingestfile, prefix=prefix, transform=transform, threshold=threshold, channel = channel, region=region)
            print("Initial sort....")
            h = sort_strains_into_histograms(data, bin_endpoints=bin_endpoints, media=media, od=od)
            print("Initial sort done.")
            cts[circuit] = h
        cts.update({'bin_centers': bin_centers})
        json.dump(cts, open(dumpfile, "w"))
        cts.pop('bin_centers')
    cts, bin_centers = format_circuits(cts,bin_centers)
    output = {}
    for circuit in circuits:
        print("Processing results for {}....".format(circuit))
        results = get_results(cts[circuit], bin_centers, num_choices=num_choices_res)
        print("Processing results done.")
        print("Building null model for {}....".format(circuit))
        null_model = build_null_model(cts[circuit], bin_centers, num_choices=num_choices_null)
        print("Building null model done.")
        output[circuit] = {'results':{str(k) : r for k,r in results.items()},'null_model':{str(k) : n for k,n in null_model.items()}}
    if not savefile:
        parts = dumpfile.split('.')
        savefile = parts[0]+"_output.json"
    json.dump(output, open(savefile, "w"))
    print("Output saved to {}".format(savefile))
    return output


def format_circuits(circuitstr,bin_centers):
    circuits = {}
    for circuit, dat in circuitstr.items():
        dats = {}
        for metadata, vals in dat.items():
            metadata = ast.literal_eval(metadata)
            vals = {k: [np.asarray(hist) for hist in vals[k]] for k in vals}
            dats[metadata] = vals
        circuits[circuit] = dats
    bin_centers = np.asarray(bin_centers)
    return circuits, bin_centers


def get_results(circuits, bin_centers, num_choices=25):
    new_circuits = {}
    for metadata, vals in circuits.items():
        md = []
        for m in metadata:
            if m[0] not in ["rep", "od"]:
                md.append(m)
        md = tuple(md)
        if md not in new_circuits:
            new_circuits[md] = vals
        else:
            temp = dict(new_circuits[md])
            new_circuits[md] = {k: temp.get(k) + vals.get(k) for k in vals.keys()}
    all_scores = {md: {'truthtable_incorrect': [], 'truthtable_correct': []} for md in new_circuits}
    for md, ip in new_circuits.items():
        for m in md:
            if m[0] == "circuit":
                desiredtt = desired_truth_tables[m[1]]
                break
        truthtable_incorrect = []
        truthtable_correct = []
        for _ in range(num_choices):
            choice = [np.asarray(random.choice(ip[input_states[k]])) for k in range(4)]
            d = dict(zip(input_states, choice))
            scores = rank_noncst_tables(d, bin_centers)
            if scores[0][1] == desiredtt:
                truthtable_correct.append(scores[1][0] - scores[0][0])
            else:
                truthtable_incorrect.append(scores[1][0] - scores[0][0])
        all_scores[md]['truthtable_incorrect'].extend(truthtable_incorrect)
        all_scores[md]['truthtable_correct'].extend(truthtable_correct)
    return all_scores


def build_null_model(circuits, bin_centers, num_choices=50):
    random_circuits = {}
    for metadata, vals in circuits.items():
        md = []
        for m in metadata:
            if m[0] not in ["rep", "od"]:
                md.append(m)
        md = tuple(md)
        all_inputs = [a for k in vals for a in vals[k]]
        if md not in random_circuits:
            random_circuits[md] = all_inputs
        else:
            random_circuits[md].extend(all_inputs)
    all_scores = {md: {'truthtable_incorrect': [], 'truthtable_correct': []} for md in random_circuits}
    for md, ip in random_circuits.items():
        for m in md:
            if m[0] == "circuit":
                desiredtt = desired_truth_tables[m[1]]
                break
        truthtable_incorrect = []
        truthtable_correct = []
        for _ in range(num_choices):
            choice = [np.asarray(random.choice(ip)) for _ in range(4)]
            d = dict(zip(input_states, choice))
            scores = rank_noncst_tables(d, bin_centers)
            if scores[0][1] == desiredtt:
                truthtable_correct.append(scores[1][0] - scores[0][0])
            else:
                truthtable_incorrect.append(scores[1][0] - scores[0][0])
        all_scores[md]['truthtable_incorrect'].extend(truthtable_incorrect)
        all_scores[md]['truthtable_correct'].extend(truthtable_correct)
    return all_scores
