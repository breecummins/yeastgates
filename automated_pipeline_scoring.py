# The MIT License (MIT)
#
# Copyright (c) 2018 Bree Cummins
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import FlowCytometryTools as FCT
import json, os
from synthetic_circuit_performance import *
from rank_order_truth_tables import *
import numpy as np
import pandas as pd
import ast, random


input_states = ["00", "01", "10", "11"]


def get_channels(lab_name):
    # key to lab specific channels
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


def transform_data(fname,transform,channels,threshold,channel,region):
    #Transform and gate using flowcytometrytools
    sample = FCT.FCMeasurement(ID="temp", datafile=fname)
    if transform == 'hlog':
        # see FlowCytometryTools documentation
        sample = sample.transform(transform, channels=[c for _, c in channels.items()], b=100)
    elif transform == 'tlog':
        sample = sample.transform(transform, channels=[c for _, c in channels.items()], th=2)
    elif transform is not None:
        raise ValueError("Transform {} not recognized.".format(transform))
    if threshold:
        try:
            gate = FCT.ThresholdGate(threshold, [channels[channel]], region=region)
            sample = sample.gate(gate)
        except:
            print(fname)
            print(sample.channel_names)
            raise
    return sample


def get_data_tx(circuit,ingest_file="transcriptic_april_fcsfiles_dan.csv", prefix="~/sd2e-community", transform='hlog',threshold=4000,channel="FSC",region="above"):
    # Parses transcriptic_april_fcsfiles_dan.csv exactly.
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
                    # wild type has no input state
                    continue
                rep = row["replicate"]
                od = row["od"]
                metadata = {'media': media, 'circuit': circuit, 'od': od, 'input_state': input_state, 'rep': rep}
                channels.pop("Sytox")
                sample = transform_data(fname, transform, channels, threshold,channel,region)
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


def bin_data(pts, bin_endpoints):
    hist = [0] * (len(bin_endpoints) + 1)
    for p in list(pts):
        ind = sorted(list(bin_endpoints) + [p]).index(p)
        hist[ind] += 1
    return hist


def get_bin_centers(bin_endpoints):
    return [a + (b - a) / 2 for (a, b) in zip([0] + list(bin_endpoints), list(bin_endpoints) + [2 * bin_endpoints[-1] - bin_endpoints[-2]])]


def get_log_values(values):
    log_vals = []
    for v in values:
        if v > 0:
            log_vals.append(np.log10(v))
        else:
            log_vals.append(-1)
    return np.asarray(log_vals)


def sort_strains_into_histograms(data, bin_endpoints):
    # Data transformation and reorganization
    circuits = {}
    for key, samp_list in data.items():
        for sample, channels, metadata in samp_list:
            ip = metadata["input_state"]
            metadata.pop("input_state")
            metadata=str(metadata)
            pts = get_log_values(sample.data[channels["GFP"]].values.transpose())
            hist = np.asarray(bin_data(pts, bin_endpoints))
            if metadata not in circuits:
                circuits[metadata] = {"00": [], "01": [], "10": [], "11": []}
            circuits[metadata][ip].append(hist)
    return circuits


def get_results(hists, bin_centers, num_choices):
    # Record separation scores and whether they are associated to the desired truth table or not.
    new_circuits = {}
    for metadata, vals in hists.items():
        metadata = ast.literal_eval(metadata)
        try:
            metadata.pop("rep")
        except:
            pass
        try:
            metadata.pop("od")
        except:
            pass
        md = tuple(metadata.items())
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
            try:
                if scores[0][1] == desiredtt:
                    truthtable_correct.append(scores[1][0] - scores[0][0])
                else:
                    truthtable_incorrect.append(scores[1][0] - scores[0][0])
            except:
                print(md)
                raise
        all_scores[md]['truthtable_incorrect'].extend(truthtable_incorrect)
        all_scores[md]['truthtable_correct'].extend(truthtable_correct)
    return all_scores


def main_tx(circuit,ingest_file="transcriptic_april_fcsfiles_dan.csv",bin_endpoints=[np.log10(r) for r in range(250, 10250, 250)], num_choices=250):
    '''
    This function works only for files in the format transcriptic_april_fcsfiles_dan.csv. There are also multiple
    default arguments in this script that came from looking at data.

    :param circuit: One of "AND", "NAND", "NOR", "OR", "XOR", "XNOR"
    :param ingest_file: see default
    :param bin_endpoints: the default is based on parameters chosen in the rest of the file
    :param num_choices: how many circuits to measure based on pooling replicates and optical densities for a fixed
    circuit, media condition, and experiment. The circuits are constructed by randomly picking each strain from the
    pooled data (pick a 00 from all the 00 strains, pick 01 from all the 01 strains, etc)
    :return: Separation scores and whether they are associated to the desired truth table are saved to a file.
    '''
    bin_centers = np.asarray(get_bin_centers(bin_endpoints))
    print("Getting data for circuit {}....".format(circuit))
    data = get_data_tx(circuit,ingest_file=ingest_file)
    print("Initial sort for circuit {}....".format(circuit))
    h = sort_strains_into_histograms(data, bin_endpoints)
    print("Initial sort done.")
    print("Processing results for {}....".format(circuit))
    results = get_results(h, bin_centers, num_choices=num_choices)
    print("Processing results done.")
    savefile = "temp_output_{}.json".format(circuit)
    json.dump({str(k) : r for k,r in results.items()}, open(savefile, "w"))
    print("Output saved to {}".format(savefile))
