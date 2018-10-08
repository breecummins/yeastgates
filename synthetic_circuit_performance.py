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


from rank_order_truth_tables import rank_noncst_tables
from math import isnan
import numpy as np
import ast,time,sys,json
import matplotlib as mpl
import matplotlib.pyplot as plt


desired_truth_tables = {
    "NOR" : {'00' : 1, '01' : 0, '10' : 0, '11' : 0},
    "OR"  : {'00' : 0, '01' : 1, '10' : 1, '11' : 1},
    "AND" : {'00' : 0, '01' : 0, '10' : 0, '11' : 1},
    "NAND": {'00' : 1, '01' : 1, '10' : 1, '11' : 0},
    "XOR" : {'00' : 0, '01' : 1, '10' : 1, '11' : 0},
    "XNOR": {'00' : 1, '01' : 0, '10' : 0, '11' : 1}
                       }


def getfilenames(linkslist):
    for l in linkslist:
        lp = l.split("tree/")
        yield "../../../" + lp[1] + "/output/output.csv"


def getcircuit(name_str):
    # order is important
    sorted_circuits = ["XNOR", "NOR", "XOR", "OR", "NAND", "AND"]
    for c in sorted_circuits:
        if c in name_str:
            return c
    return None

def getinputstate(name_str):
    inputs = ["00", "10", "01", "11"]
    for i in inputs:
        if i in name_str:
            return i
    return None

def printfile(fname):
    fn = fname.split('/')
    N = len(fn) // 2 
    print('/'.join(fn[:N])+'/\n'+'/'.join(fn[N:]))

def parseoutput(f,underflow_skip=2,truncate_ind=1):
    '''
    This parses a specific data product from TASBE which is likely obsolete.

    :param f: file object
    :param underflow_skip: skip this number of bins at the beginning
    :param truncate_ind: skip this number of bins at the end
    :return: dictionary of specific metadata each associated to a histogram and the representative bin values
    '''
    out = {}
    f0 = f.readline()
    l = f0.split('geo_mean,')
    bin_vals = np.array([float(L) for L in l[1].split(',')][underflow_skip:-truncate_ind]) #skip underflow channel, truncate high bins
    for l in f:
        u=l.split('}",')
        nums = u[1].split(',')
        c = ast.literal_eval(u[0][1:] + '}')
        g = float(nums[1])
        d = np.array([float(N) for N in nums[2+underflow_skip:-truncate_ind]]) #skip underflow channel, truncate high bins
        try:
            strain = c['design_name']
            circuit = getcircuit(strain)
            if circuit:
                inputstate = getinputstate(strain)
                rep = "Replicate " + str(c["replicate"]) + ", target od = " + str(c["target_od"])
                media = c["media_name"] 
                try:
                    out[rep].update({ (circuit,inputstate,media,g) : d })
                except:
                    out.update({rep : { (circuit,inputstate,media,g) : d }})
        except:
            pass
    return out, bin_vals


def make_hists(data,bin_vals,desiredtt,rep):
    plt.figure()
    plt.title("{}".format(rep))
    for k,hist in data.items():
        if desiredtt[k] == 1:
            c = 'blue'
        else:
            c = 'red'
        H1 = hist / float(sum(hist))
        mpl.rc('xtick', labelsize=14) 
        mpl.rc('ytick', labelsize=14)         
        plt.hist(bin_vals,len(bin_vals),weights=H1,alpha=0.5,color=c)
    plt.show()


def do_all(fname, good_sep_param, print_reps=True,print_summary=True,make_figs=True,gates=["OR","NOR","AND","NAND","XOR","XNOR"],top_tables=7):
    '''
    This function is obsolete, both because of parsing and because now I am mixing across replicates.

    :param fname: the file to parse
    :param good_sep_param: an arbitrary choice of a good separation score
    :param print_reps: print info as we go
    :param print_summary: print info as we go
    :param make_figs: make figures as we go
    :param gates: any sublist of the default argument
    :param top_tables: integer between 0 and 7: the number of truth tables to print results for
    :return:
    '''
    starttime=time.time()
    with open(fname,'r') as f:
        out, bin_vals = parseoutput(f)
        circuits = list(set([tup[0] for rep in out for tup in out[rep]]))
        if len(circuits) > 1:
            printfile(fname)
            print(circuits)
            raise ValueError("Input file has multiple circuits. Inspect for parsing.")
        elif circuits[0] in gates:
            printfile(fname)
            correct = 0
            goodsets = 0
            good_sep = 0
            savescores = []
            for rep, hist in out.items():
                keys = [ tup for tup in out[rep] ]
                circuit =  keys[0][0]
                if circuit in gates:
                    if print_reps:
                        print((circuit,keys[0][2]))
                    # Don't process if there are nan's
                    if any([isnan(k[3]) for k in keys]):
                        print("Missing data for conditions {}.".format(" ".join(k[1] for k in keys if isnan(k[3]))))
                        continue
                    desiredtt = desired_truth_tables[circuit]
                    data = { k[1] : out[rep][k] for k in keys }
                    scores = rank_noncst_tables( data, bin_vals )
                    goodsets += 1 
                    if scores[0][1] == desiredtt:
                        correct += 1
                        savescores.append((scores[0][0], scores[1][0] - scores[0][0],True))
                        if scores[1][0] - scores[0][0] > good_sep_param:
                            good_sep += 1
                    else:
                        savescores.append((scores[0][0], scores[1][0] - scores[0][0], False))
                    if make_figs:
                        make_hists(data,bin_vals,desiredtt,rep)
                        time.sleep(0.1)
                    if print_reps:
                        print(rep)
                        print_scores(scores[:top_tables])
            if print_summary:
                print((circuit,keys[0][2]))
                if goodsets > 0:
                    print("{} correct / {} replicates = {:.2f}%, {} with good separation / {} replicates = {:.2f}%".format(correct,goodsets,float(correct)/goodsets*100,good_sep,goodsets,float(good_sep)/goodsets*100))
                else:
                    print("No good replicates.")
                print(time.time()-starttime)
                print("\n")
                print("---------------------------------------------------------------------------------")
                print("---------------------------------------------------------------------------------")
                print("---------------------------------------------------------------------------------")
                print("\n")
            return savescores,circuit,keys[0][2]
