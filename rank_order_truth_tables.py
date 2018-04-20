import pyemd, itertools
import numpy as np


def make_bin_dist(bin_vals):
    '''
    Make a distance matrix for the bin values of a histogram.

    :param bin_vals: list or 1-D numpy array of representative points (usually midpts) of the bins defining a collection
    of histograms. For yeast gates, these values are generally going to be exponents of fluorescence.
    :return: a numpy array of the pairwise distances between points
    '''
    return np.array([[np.abs(m - n) for m in bin_vals] for n in bin_vals])


def emdist(h1, h2, bin_dist):
    '''
    Calculate earth mover's distance between 2 histograms.

    :param h1: a 1-D numpy array representing a histogram with the bin values used to make bin_hist
    :param h2: a 1-D numpy array representing a histogram with the bin values used to make bin_hist
    :param bin_dist: output from make_bin_dist(bin_vals)
    :return: a scalar value that is the earth mover's distance between normalized h1 and h2
    '''
    return pyemd.emd(h1/float(sum(h1)), h2/float(sum(h2)), bin_dist)


def similarity(dist):
    '''
    Transform a distance into a similarity score.

    :param dist: a non-negative scalar
    :return: a scalar between 0 and 1
    '''
    return np.exp(-dist**2 / 2)


def bestcontrol(inputstates, data, hmin, hmax, bin_dist):
    '''
    Calculate which control distribution is closer to the group of histograms determined by the inputs in inputstates.

    :param inputstates: a subset of the keys of data
    :param data: a length 2^n dictionary mapping an input state to a 1-D numpy array histogram, all defined by the same bin values
    :param hmin: the minimal control histogram, 1-D numpy array
    :param hmax: the maximal control histogram, 1-D numpy array
    :param bin_dist: output from make_bin_dist(bin_vals), where bin_vals are the bin values defining the histograms in data
    :return: 'min' or 'max' denoting the closest control histogram
    '''
    vmin, vmax = 0, 0
    for i in inputstates:
        vmin += emdist(data[i], hmin, bin_dist)
        vmax += emdist(data[i], hmax, bin_dist)
    return 'min' if vmin <= vmax else 'max'


def calculate_partitions(N):
    m = int(N/2)
    partitions = []
    for k in range(m, N+1):
        combs = [c for c in itertools.combinations(range(N), k)]
        if k == m:
            # there is double-counting of partitions at the median value
            partitions.extend(combs[:len(combs)/2])
        else:
            partitions.extend(combs)
    return partitions


def make_truth_tables(N, partitions):
    truthtables = []
    for p in partitions:
        truthtables.append([0 if a in p else 1 for a in range(N)])
        truthtables.append([1 if a in p else 0 for a in range(N)])
    return truthtables


def rank_by_median_similarity(data):
    inputstates = [k for k in data]
    partitions = calculate_partitions(len(inputstates))
    truthtables = make_truth_tables(len(inputstates), partitions)
