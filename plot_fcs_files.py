import pylab

def plot_2D_scatter(data, xchannel="GFP", ychannel="FSC", xlim=[-0.01, 7], ylim=[-0.01, 7], pool=False, savefile=False):
    '''

    :param data: the output of get_data() of process_fcs_files
    :param xchannel: one of ["GFP","FSC","Sytox"]
    :param ychannel: one of ["GFP","FSC","Sytox"]
    :param xlim: plot limits
    :param ylim: plot limits
    :param pool: pool data across media conditions
    :return:
    '''
    for t, d in data.items():
        media = set([])
        for sample, channels, metadata in d:
            # alter plotting for log10
            sample.plot([channels[xchannel], channels[ychannel]], cmap='hot')
            media.add(metadata['media'])
        if not pool:
            title = media[0] if len(media) == 1 else t
            pylab.title(title)
            pylab.xlabel(xchannel)
            pylab.ylabel(ychannel)
            pylab.xlim(xlim)
            pylab.ylim(ylim)
            if savefile:
                pylab.savefig(savefile)
            pylab.show()
    if pool:
        title = t
        pylab.title(title)
        pylab.xlabel(xchannel)
        pylab.ylabel(ychannel)
        pylab.xlim(xlim)
        pylab.ylim(ylim)
        if savefile:
            pylab.savefig(savefile)
        pylab.show()


def plot_hist(data,xchannel="FSC",xlim=[6500,12000],ylim=[0,1400],bins=200,color="blue",pool=False,savefile=False):
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
        media = set([])
        for sample,channels,metadata in d:
            sample.plot([channels[xchannel]], bins=bins,color=color)
            media.add(metadata['media'])
        if not pool:
            title = media[0] if len(media) == 1 else t
            pylab.title(title)
            pylab.xlabel(xchannel)
            pylab.ylabel(ychannel)
            pylab.xlim(xlim)
            pylab.ylim(ylim)
            if savefile:
                pylab.savefig(savefile)
            pylab.show()
    if pool:
        title = t
        pylab.title(title)
        pylab.xlabel(xchannel)
        pylab.ylabel(ychannel)
        pylab.xlim(xlim)
        pylab.ylim(ylim)
        if savefile:
            pylab.savefig(savefile)
        pylab.show()


def plot_hist_from_point_cloud(ptclouds,xlabel="",title="",bin_endpts=None,xlim=None,ylim=None,colors=["blue"],bins=20,savefile=False):
    # ptclouds is list of np.arrays or lists, each for a different input condition
    # all will be plotted in the same figure
    for ptcloud,color in zip(ptclouds,colors):
        if bin_endpts is not None:
            pylab.hist(ptcloud,bins=bin_endpts,color=color,alpha=.5)
        else:
            pylab.hist(ptcloud,bins=bins,color=color,alpha=.5)
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


def plot_bar(x,heights,xlabel="",ylabel="counts",title="",xlim=[0,0.1],ylim=[0,50],colors=["blue"],savefile=False):
    for h,c in zip(heights,colors):
        pylab.hist(x, len(x), weights=h/sum(h), alpha=0.5, color=c)
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    if xlim:
        pylab.xlim(xlim)
    if ylim:
        pylab.ylim(ylim)
    if savefile:
        pylab.savefig(savefile)
    pylab.show()
    
    
def sep_hists(results,title,scale=10,bins=50):
    # results is the output of process_fcs_files.get_results 
    ptcloud = [scale*r for r in results[list(results.keys())[0]]["truthtable_incorrect"]]
    try:
        M=max(ptcloud)
    except:
        M=1    
    bin_endpts = np.linspace(0,M,bins)
    plotff.plot_hist_from_point_cloud([ptcloud],xlim=[0,M],ylim=[0,len(ptcloud)],title=title,xlabel="incorrect separation scores * {}".format(scale),bin_endpts=bin_endpts,colors=["red"])

    ptcloud = [scale*r for r in results[list(results.keys())[0]]["truthtable_correct"]]
    try:
        M=max(ptcloud)
    except:
        M=1    
    bin_endpts = np.linspace(0,M,bins)
    plotff.plot_hist_from_point_cloud([ptcloud],xlim=[0,M],ylim=[0,len(ptcloud)],title=title,xlabel="correct separation scores * {}".format(scale),bin_endpts=bin_endpts,colors=["green"])



def plot_summaries(results,ylim=(0,0.04)):
    # results is the output of process_fcs_files.get_results 
    circuits = ["OR","NOR","AND","NAND","XOR","XNOR"]
    input_states = ["00","01","10","11"]
    SC = {}
    Sorb = {}
    Eth = {}
    YEPD = {}
    media_scores = {"SC":SC,"Sorb":Sorb,"Eth":Eth,"YEPD":YEPD}
    for m,scores in results.items():
        circuit = pff.getcircuit(m[1][1])
        media = m[0][1]
        if media == "Yeast_Extract_Peptone_Adenine_Dextrose" or media == "culture_media_4":
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
        plt.figure()
        plt.title(m)
        ax = plt.gca()
        ax.set_xticks(range(6))
        ax.set_xticklabels(circuits, fontsize=18)
        plt.ylim(ylim)
        for c,cor in cts.items():
            ind = circuits.index(c)
            seps = cor["truthtable_correct"]
            inds = [ind]*len(seps)
            col = "g"
            plt.plot(inds,seps,color=col,marker="o")
            seps = cor["truthtable_incorrect"]
            inds = [ind]*len(seps)
            col = "r"
            plt.plot(inds,seps,color=col,marker="o")
        plt.show()
        


