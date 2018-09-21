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


def plot_hist_from_point_cloud(ptclouds,xlabel="",title="",bin_endpts=None,xlim=[0,0.1],ylim=[0,50],colors=["blue"],bins=20,savefile=False):
    for ptcloud,color in zip(ptclouds,colors):
        if bin_endpts:
            pylab.hist(ptcloud,bins=bin_endpts,color=color)
        else:
            pylab.hist(ptcloud,bins=bins,color=color)
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel("count")
    pylab.xlim(xlim)
    pylab.ylim(ylim)
    if savefile:
        pylab.savefig(savefile)
    pylab.show()
