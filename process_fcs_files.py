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
from rank_order_truth_tables import rank_noncst_tables
from synthetic_circuit_performance import getcircuit
import numpy as np
import itertools
import pandas as pd
import ast, random
from pprint import pprint

input_states = ["00", "01", "10", "11"]

desired_truth_tables = {
    "NOR" : {'00' : 1, '01' : 0, '10' : 0, '11' : 0},
    "OR"  : {'00' : 0, '01' : 1, '10' : 1, '11' : 1},
    "AND" : {'00' : 0, '01' : 0, '10' : 0, '11' : 1},
    "NAND": {'00' : 1, '01' : 1, '10' : 1, '11' : 0},
    "XOR" : {'00' : 0, '01' : 1, '10' : 1, '11' : 0},
    "XNOR": {'00' : 1, '01' : 0, '10' : 0, '11' : 1}
                       }


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


def transform_data(fname,transform,channels,threshold,channel,region):
    sample = FCT.FCMeasurement(ID="temp", datafile=fname)
    if transform == 'hlog':
        # see FlowCytometryTools documentation
        sample = sample.transform(transform, channels=[c for _, c in channels.items()], b=100)
    elif transform == 'tlog':
        sample = sample.transform(transform, channels=[c for _, c in channels.items()], th=2)
    elif transform is not None:
        raise ValueError("Transform {} not recognized.".format(transform))
    if threshold:
        gate = FCT.ThresholdGate(threshold, [channels[channel]], region=region)
        sample = sample.gate(gate)
    return sample

def get_data(circuit,ingest_file="matches_biofab_all_circuits_all_media.json",transform=None,threshold=None,channel=None,region=None):
    matches = json.load(open(ingest_file))
    data = {}
    count = 0
    for m in matches:
        try:
            c = getcircuit(m['strain_circuit'])
        except:
            continue
        if c == circuit:
            lab_name = m['lab']
            channels = get_channels(lab_name)
            fname = m['jupyter_path']
            media = m['sample_contents'][0]['name']['label']
            input_state = m['strain_input_state']
            try:
                rep = m['replicate']
            except:
                rep = None
            try:
                od = m['inoculation_density']['value']
            except:
                od = None
            metadata = {'media': media, 'circuit': circuit, 'od': od, 'input_state': input_state, 'rep': rep}
            sample = transform_data(fname, transform, channels, threshold,channel,region)
            d = fname.split('/')[-2]
            if d in data:
                data[d].append((sample, channels, metadata))
            else:
                data[d] = [(sample, channels, metadata)]
            count += 1
            if not count % 50:
                print("{} files found.".format(count))
    print("{} files in total.".format(count))
    return data


def get_log_values(values):
    log_vals = []
    for v in values:
        if v > 0:
            log_vals.append(np.log10(v))
        else:
            log_vals.append(-1)
    return np.asarray(log_vals)


def bin_data(pts, bin_endpoints):
    hist = [0] * (len(bin_endpoints) + 1)
    for p in list(pts):
        ind = sorted(list(bin_endpoints) + [p]).index(p)                
        hist[ind] += 1
    return hist


def get_bin_centers(bin_endpoints):
    return [a + (b - a) / 2 for (a, b) in zip([0] + list(bin_endpoints), list(bin_endpoints) + [2 * bin_endpoints[-1] - bin_endpoints[-2]])]


def sort_strains_into_histograms(data, bin_endpoints, media=None, od=None):
    circuits = {}
    for key, samp_list in data.items():
        for sample, channels, metadata in samp_list:
            if (not media or metadata["media"] == media) and (not od or metadata['od'] == od):
                print(metadata)
                ip = metadata["input_state"]
                md = metadata.copy()
                md.pop("input_state")
                md = tuple(md.items())
                pts = get_log_values(sample.data[channels["GFP"]].values.transpose())
                hist = bin_data(pts, bin_endpoints)
                if md not in circuits:
                    circuits[md] = {"00": [], "01": [], "10": [], "11": []}
                circuits[md][ip].append(np.asarray(hist))
    return circuits


def get_results(circuits, bin_endpoints, num_choices=25):
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
    bin_centers = np.asarray(get_bin_centers(bin_endpoints))
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


def build_null_model(circuits, bin_endpoints, num_choices=50):
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
    bin_centers = get_bin_centers(bin_endpoints)
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
