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


def emdist(h1, h2, bin_vals):
    '''
    Calculate earth mover's distance between 2 histograms.

    :param h1: a 1-D numpy array representing a histogram with the bin values used to make bin_dist
    :param h2: a 1-D numpy array representing a histogram with the bin values used to make bin_dist
    :param bin_vals: representative points in the bins defining the histograms in data
    :return: a scalar value that is the earth mover's distance between normalized h1 and h2
    '''
    return pyemd.emd(h1/float(sum(h1)), h2/float(sum(h2)), make_bin_dist(bin_vals))


def similarity(dist):
    '''
    Transform a distance into a similarity score.

    :param dist: a non-negative scalar
    :return: a scalar between 0 and 1
    '''
    return np.exp(-dist**2 / 2)


def make_graph(data,bin_vals):
    '''
    Make a dictionary of pairwise similarity values that represents a weighted graph.

    :param data: a length 2^k dictionary mapping an input state to a 1-D numpy array histogram, all defined by the same bin values
    :param bin_vals: representative points in the bins defining the histograms in data
    :return: dictionary of real values between 0 and 1 keyed by pairs of input states, with no repeated pairs (upper triangular part of similarity matrix, or similarly, the list of weighted edges in an undirected graph)
    '''
    inputstates = [k for k in data]
    graph = {}
    for t,i in enumerate(inputstates[-1]):
        for j in inputstates[t+1:]:
            graph.update( {(i,j) : similarity(emdist(data[i],data[j],bin_vals))} )
    return graph


def mean_hist(hist,bin_vals):
    '''
    Return the mean of a histogram.

    :param hist: an iterable containing the counts in each bin
    :param bin_vals: an iterable containing a representative value for each bin
    :return: mean value of the histogram
    '''
    return sum([x*h for (x,h) in zip(bin_vals,hist)]) / float(sum(hist))


def calculate_partitions(inputstates):
    '''
    Calculate all ways to partition a list into 2 nonempty sets, where the list has an even number of elements.

    :param inputstates: a list with an even number of elements, N
    :return: a length 2^(N-1) - 1 list of lists, where each sublist represents half of a partition (the other half is deduced by inputstates \ sublist)
    '''
    N = len(inputstates)
    partitions = []
    for k in range(int(N/2)+1):
        combs = [c for c in itertools.combinations(inputstates, k)]
        partitions.extend(combs)
    # there is double-counting of partitions at the median value
    partitions = partitions[:-int(len(combs)/2)]
    return partitions


def make_truth_table(partition,inputstates,b):
    '''
    Return the truth table for a partition.

    :param partition: an iterable containing a subset of inputstates
    :param inputstates: a list of length N
    :param b: 0 or 1
    :return: a length N list of 0's and 1's representing a truth table
    '''
    return [(i,b) if i in partition else (i, (b + 1) % 2) for i in inputstates]


def weight_to_part(p, graph):
    '''
    Calculate the weight function for the normalized cut score.

    :param p: one set in a 2-partition of inputstates
    :param graph: the similarity values for the inputstates graph
    :return: a scalar value
    '''
    return sum( [ v for (i,j),v in graph.items() if (i in p) != (j in p) ] )


def weight_to_all(p, graph):
    '''
    Calculate the weight function for the normalized cut score.

    :param p: one set in a 2-partition of inputstates
    :param graph: the similarity values for the inputstates graph
    :return: a scalar value
    '''
    return sum( [ v for (i,j),v in graph.items() if i in p ] )


def normalized_cut(partition, inputstates, graph):
    '''
    Calculate normalized cut score for the partition and inputstates \ partition.

    :param partition: a list or set containing a subset of inputstates
    :param inputstates: a list or set of input states (the keys to data)
    :param graph: a dictionary that is the output of make_graph()
    :return: a scalar value
    '''
    Wpq = weight_to_part(partition,graph)
    WpV = weight_to_all(partition,graph)
    WqV = weight_to_all(set(inputstates).difference(partition),graph)
    return Wpq * (1.0/WpV + 1.0/WqV)

def rank_nonconstant_tables(data,bin_vals):
    '''
    Return scored nonconstant truth tables.

    :param data: a length 2^k dictionary mapping an input state to a 1-D numpy array histogram, all defined by the same bin values
    :param bin_vals: representative points in the bins defining the histograms in data
    :return: dictionary of real values keyed by truth tables
    '''
    graph = make_graph(data,bin_vals)
    inputstates = [k for k in data]
    partitions = calculate_partitions(inputstates)
    scores = {}
    for p in partitions:
        hp = sum( [ h for k,h in data.items() if k in p ] )
        hq = sum( [ h for k,h in data.items() if k not in p ] )
        mp = mean_hist(hp,bin_vals)
        mq = mean_hist(hq,bin_vals)
        b = 0 if mp < mq else 1
        truthtable = make_truth_table(p,inputstates,b)
        scores.update( { truthtable : normalized_cut(p, inputstates, graph) } )
    return scores


def print_scores(scores):
    for tt,ns in scores.items():
        print(ns)
        print("\n")
        for t in tt:
            print(t)
        print("\n")

def assess_cst_truthtable():
    pass
