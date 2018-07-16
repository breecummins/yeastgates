from rank_order_truth_tables import *
import numpy as np
import ast,time,sys,json
import matplotlib as mpl
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


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

def parseoutput(f,underflow_skip=2,truncate_ind=1):    
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


def do_all(fname, good_sep_param, print_reps=True,make_figs=True,gates=["OR","NOR","AND","NAND","XOR","XNOR"],top_tables=7):
    starttime=time.time()
    with open(fname,'r') as f:
        out, bin_vals = parseoutput(f)
        circuits = list(set([tup[0] for rep in out for tup in out[rep]]))
        if len(circuits) > 1:
            print(circuits)
            raise ValueError("Input file has multiple circuits. Inspect for parsing.")
        elif circuits[0] in gates:
            print(fname)
            correct = 0
            goodsets = 0
            good_sep = 0
            for rep, hist in out.items():
                keys = [ tup for tup in out[rep] ]
                circuit =  keys[0][0]
                if circuit in gates:
                    if print_reps:
                        print((circuit,keys[0][2]))
                    desiredtt = desired_truth_tables[circuit]
                    data = { k[1] : out[rep][k] for k in keys }
                    scores = rank_noncst_tables( data, bin_vals )
                    if scores[0][0] > 0: #filter out nan's
                        goodsets += 1 
                        if scores[0][1] == desiredtt:
                            correct += 1
                            if scores[1][0] - scores[0][0] > good_sep_param:
                                good_sep += 1
                    if make_figs:
                        make_hists(data,bin_vals,desiredtt,rep)
                        time.sleep(0.1)
                    if print_reps:
                        print(rep)
                        print_scores(scores[:top_tables])
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
