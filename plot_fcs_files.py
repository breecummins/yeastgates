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


import pylab
import numpy as np
import process_fcs_files as pff
import ast, random


def plot_2D_scatter(data, media = "",titlestr="",xchannel="GFP", ychannel="FSC", xlim=None, ylim=None, savefile=None):
    '''

    :param data: the output of get_data() of process_fcs_files
    :param xchannel: one of ["GFP","FSC","Sytox"]
    :param ychannel: one of ["GFP","FSC","Sytox"]
    :param xlim: plot limits
    :param ylim: plot limits
    :return:
    '''
    for t, d in data.items():
        med_list=set([])
        for sample, channels, metadata in d:
            if (not media) or metadata['media'] == media:
                sample.plot([channels[xchannel], channels[ychannel]], cmap='hot')
                med_list.add(metadata['media'])
        title = titlestr+ " " + t + " "+ " & ".join(med_list)
        pylab.title(title)
        pylab.xlabel(xchannel)
        pylab.ylabel(ychannel)
        if xlim:
            pylab.xlim(xlim)
        if ylim:
            pylab.ylim(ylim)
        if savefile:
            pylab.savefig(savefile)
        pylab.show()
        
        
def plot_data_against_gated(circuit,mongo=False, ingest_file="transcriptic_april_fcsfiles_dan.csv", prefix="~/sd2e-community/",transform='hlog',threshold=5000,channel = "FSC",region="above", xlim=[-7000, 10000], ylim=[-7000, 10000],Sytox=False):
    
    def make_thresh(data,threshold):
        if threshold:
            gated_data = {}
            for d,l in data.items():
                new_l = []
                for (sample,channels,metadata) in l:
                    gated_sample = pff.make_threshold_gate(sample,threshold,channels[channel],region=region)
                    new_l.append((gated_sample,channels,metadata))
                gated_data[d] = new_l
            data = gated_data
        return data
 
    if mongo:
        data = pff.get_data_mongo(circuit,ingest_file, prefix, transform)
    else:
        data = pff.get_data_tx(circuit,ingest_file, prefix, transform)
      
    plot_2D_scatter(data, xchannel="GFP", ychannel="FSC", xlim=xlim, ylim=ylim, pool=False, savefile=False)
    if Sytox:
        plot_2D_scatter(data, xchannel="GFP", ychannel="Sytox", xlim=xlim, ylim=ylim, pool=False, savefile=False)
    
    if threshold:
        gdata = make_thresh(data,threshold)
        plot_2D_scatter(gdata, titlestr = "gated", xchannel="GFP", ychannel="FSC", xlim=xlim, ylim=ylim, pool=False, savefile=False)
        if Sytox:
            plot_2D_scatter(gdata, titlestr = "gated", xchannel="GFP", ychannel="Sytox", xlim=xlim, ylim=ylim, pool=False, savefile=False)
        data = gdata
#     return data

    
def plot_hist(data,xchannel="FSC",xlim=None,ylim=None,bins=200,color="blue",savefile=False):
    '''

    :param data: the output of get_data() of process_fcs_files
    :param xchannel: one of ["GFP","FSC","Sytox"]
    :param ychannel: one of ["GFP","FSC","Sytox"]
    :param xlim: plot limits
    :param ylim: plot limits
    :param pool: pool data across media conditions
    :return:
    '''
    for t,d in data.items():
        med_list=set([])
        for sample, channels, metadata in d:
            if (not media) or metadata['media'] == media:
                sample.plot([channels[xchannel]], bins=bins,color=color)
                med_list.add(metadata['media'])
        title = titlestr+ " " + t + " "+ " & ".join(med_list)
        pylab.title(title)
        pylab.xlabel(xchannel)
        pylab.ylabel(ychannel)
        if xlim:
            pylab.xlim(xlim)
        if ylim:
            pylab.ylim(ylim)
        if savefile:
            pylab.savefig(savefile)
        pylab.show()
        

def plot_hist_from_point_cloud(ptclouds,xlabel="",title="",bin_endpts=None,xlim=None,ylim=None,colors=["blue"],bins=20,savefile=False,normed=0):
    # ptclouds is list of np.arrays or lists, each for a different input condition
    # all will be plotted in the same figure
    for ptcloud,color in zip(ptclouds,colors):
        if bin_endpts is not None:
            pylab.hist(ptcloud,bins=bin_endpts,color=color,alpha=.5,normed=normed)
        else:
            pylab.hist(ptcloud,bins=bins,color=color,alpha=.5,normed=normed)
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel("count")
    if xlim:
        pylab.xlim(xlim)
    if ylim:
        pylab.ylim(ylim)
    if savefile:
        pylab.savefig(savefile)
    pylab.show()


def make_hists(data,bin_vals,desiredtt,rep,circuit):
    plt.figure()
    plt.title("{}, {}".format(circuit,rep))
    plt.xlabel("GFP")
    plt.ylabel("normalized counts")
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


def plot_example_circuit(circuit,hists,bin_endpoints,xlabel="",ylabel="normalized counts",num_examples=1,xlim=None,ylim=None,savefile=False):
    desiredtt = pff.desired_truth_tables[circuit]
    colors = ["b"  if desiredtt[ip] == 0  else "r" for ip in pff.input_states]
    bin_centers = pff.get_bin_centers(bin_endpoints)
    combined_hists = pff.combine_over_rep_od(hists)
    for md,ipdict in combined_hists.items():
        for _ in range(num_examples):
            example = [np.asarray(random.choice(ipdict[pff.input_states[k]])) for k in range(4)]
            pylab.figure()
            pylab.title(str(md))
            pylab.xlabel(xlabel)
            pylab.ylabel(ylabel)
            for h,c in zip(example,colors):
                pylab.hist(bin_centers, len(bin_centers), weights=h/sum(h), alpha=0.5, color=c)
            if xlim:
                pylab.xlim(xlim)
            if ylim:
                pylab.ylim(ylim)
            if savefile:
                pylab.savefig(savefile)
            pylab.show()
    
    
def sep_hists(results,title,scale=10,bins=20,M=None):
    # results is the output of process_fcs_files.get_results 
    for key,res in results.items():
        ptcloud = [scale*r for r in res["truthtable_incorrect"]]
        if not M:
            try:
                M=max(ptcloud)
            except:
                M=1    
        bin_endpts = np.linspace(0,M,bins)[1:]
        hist = pff.bin_data(ptcloud,bin_endpts)
        bin_centers = pff.get_bin_centers(bin_endpts)
        make_hist(hist,bin_centers,xlim=[0,M],ylim=[0,1],title=title+" "+str(key),xlabel="incorrect separation scores * {}".format(scale),color="red")

        ptcloud = [scale*r for r in res["truthtable_correct"]]
        if not M:
            try:
                M=max(ptcloud)
            except:
                M=1    
        bin_endpts = np.linspace(0,M,bins)
        hist = pff.bin_data(ptcloud,bin_endpts)
        bin_centers = pff.get_bin_centers(bin_endpts)
        make_hist(hist,bin_centers,xlim=[0,M],ylim=[0,1],title=title+" "+str(key),xlabel="correct separation scores * {}".format(scale),color="green")


def make_hist(hist,bin_vals,normed=True,xlim=None,ylim=None,title="",xlabel="",ylabel="",color="blue"):
    pylab.figure()
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    if normed:
        hist = np.asarray(hist) / float(sum(hist))
    if xlim:
        pylab.xlim(xlim)
    if ylim:
        pylab.ylim(ylim)
#     mpl.rc('xtick', labelsize=14) 
#     mpl.rc('ytick', labelsize=14) 
    pylab.hist(bin_vals,len(bin_vals),weights=hist,alpha=0.5,color=color,)
    pylab.show()


def plot_summaries(results,ylim=(0,0.04)):
    # results is the output of process_fcs_files.get_results 
    circuits = ["OR","NOR","AND","NAND","XOR","XNOR"]
    SC = {}
    Sorb = {}
    Eth = {}
    YEPD = {}
    media_scores = {"SC":SC,"Sorb":Sorb,"Eth":Eth,"YEPD":YEPD}
    for m,scores in results.items():
        m = ast.literal_eval(m)
        circuit = pff.getcircuit(m[1][1])
        media = m[0][1]
        if media == "Yeast_Extract_Peptone_Adenine_Dextrose" or media == "culture_media_4" or media == "YPAD":
            media = "YEPD"
        elif media == "Synthetic_Complete" or media == "culture_media_5" or media == "culture_media_1":
            media = "SC"
        elif "thanol" in media or media == "culture_media_2":
            media = "Eth"
        elif "orbitol" in media or media == "culture_media_3":
            media = "Sorb"
        else:
            raise ValueError("Media {} is not recognized.".format(media))
        media_scores[media].update({circuit : scores})
    for m,cts in media_scores.items():
        pylab.figure()
        pylab.title(m)
        ax = pylab.gca()
        ax.set_xticks(range(6))
        ax.set_xticklabels(circuits, fontsize=18)
        pylab.ylim(ylim)
        for c,cor in cts.items():
            ind = circuits.index(c)
            seps = cor["truthtable_correct"]
            inds = [ind]*len(seps)
            col = "g"
            pylab.plot(inds,seps,color=col,marker="o",linestyle="")
            seps = cor["truthtable_incorrect"]
            inds = [ind]*len(seps)
            col = "r"
            pylab.plot(inds,seps,color=col,marker="o",linestyle="")
        pylab.show()
        

def plot_by_media(circuit,results,ylim=None):
    all_media = ["SC","Sorb","YEPD","Eth"]
    fig=pylab.figure()
    pylab.title(circuit)
    ax = pylab.gca()
    pylab.xlim((-0.5,3.5))
    ax.set_xticks(range(4))
    ax.set_xticklabels(all_media, fontsize=18)
    if ylim:
        pylab.ylim(ylim)
    pylab.ylabel("separation score")

    for m,scores in results.items():
        media = m[0][1]
        if media == "Yeast_Extract_Peptone_Adenine_Dextrose" or media == "culture_media_4" or media == "YPAD":
            media = "YEPD"
        elif media == "Synthetic_Complete" or media == "culture_media_5" or media == "culture_media_1": # media 5 mistake in April data
            media = "SC"
        elif "thanol" in media or media == "culture_media_2":
            media = "Eth"
        elif "orbitol" in media or media == "culture_media_3":
            media = "Sorb"
        else:
            raise ValueError("Media {} is not recognized.".format(media))
        ind = all_media.index(media)
        vals = scores["truthtable_incorrect"]
        pylab.plot([ind]*len(vals),vals, color = "r",marker="o",linestyle="")
        vals = scores["truthtable_correct"]
        pylab.plot([ind]*len(vals),vals, color = "g",marker="o",linestyle="")
    pylab.show()

