import pyemd, itertools
import numpy as np
from pprint import pprint


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
    return pyemd.emd(np.asarray(h1)/float(sum(h1)), np.asarray(h2)/float(sum(h2)), make_bin_dist(bin_vals))


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
    for t,i in enumerate(inputstates[:-1]):
        for j in inputstates[t+1:]:
            v = emdist(data[i],data[j],bin_vals)
            graph.update( {(i,j) : similarity(emdist(data[i],data[j],bin_vals))} )
    return graph


def calculate_partitions(inputstates):
    '''
    Calculate all ways to partition a list into 2 nonempty sets, where the list has an even number of elements.

    :param inputstates: a list with an even number of elements, N
    :return: a length 2^(N-1) - 1 list of lists, where each sublist represents half of a partition (the other half is deduced by inputstates \ sublist)
    '''
    N = len(inputstates)
    partitions = []
    for k in range(1,int(N/2)+1):
        combs = [c for c in itertools.combinations(inputstates, k)]
        partitions.extend(combs)
    # there is double-counting of partitions at the median value
    partitions = partitions[:-int(len(combs)/2)]
    return partitions


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
    return sum( [ v for (i,j),v in graph.items() if (i in p) or (j in p) ] )


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


def mean_hist(hist,bin_vals):
    '''
    Return the mean of a histogram.

    :param hist: an iterable containing the counts in each bin
    :param bin_vals: an iterable containing a representative value for each bin
    :return: mean value of the histogram
    '''
    return sum([x*h for (x,h) in zip(bin_vals,hist)]) / float(sum(hist))


def assign_0or1(partition,data,bin_vals):
    '''
    Figure out if a partition should be assigned 1 (higher valued collection of histograms) or 0 (lower valued)

    :param partition: a list or set containing a subset of the keys in data
    :param data: a length 2^k dictionary mapping an input state to a 1-D numpy array histogram, all defined by the same bin values
    :param bin_vals: representative points in the bins defining the histograms in data
    :return: 0 or 1
    '''
    hp = sum([h for k, h in data.items() if k in partition])
    hq = sum([h for k, h in data.items() if k not in partition])
    mp = mean_hist(hp, bin_vals)
    mq = mean_hist(hq, bin_vals)
    # hp = [mean_hist(h,bin_vals) for k,h in data.items() if k in partition]
    # hq = [mean_hist(h,bin_vals) for k,h in data.items() if k not in partition]
    # mp  = sum(hp) / len(hp)
    # mq  = sum(hq) / len(hq)
    return 0 if mp < mq else 1


def make_truth_table(partition,inputstates,b):
    '''
    Return the truth table for a partition.

    :param partition: an iterable containing a subset of inputstates
    :param inputstates: a list of length N
    :param b: 0 or 1
    :return: a dictionary of input states keying Boolean values
    '''
    return dict([(i,b) if i in partition else (i, (b + 1) % 2) for i in inputstates])


def rank_noncst_tables(data,bin_vals):
    '''
    Return scored non-constant truth tables.

    :param data: a length 2^k dictionary mapping an input state to a 1-D numpy array histogram, all defined by the same bin values
    :param bin_vals: representative points in the bins defining the histograms in data
    :return: list of tuples of a real valued normalized cut score with its associated truth table
    '''
    graph = make_graph(data,bin_vals)
    inputstates = [k for k in data]
    partitions = calculate_partitions(inputstates)
    scores = []
    for p in partitions:
        truthtable = make_truth_table(p,inputstates,assign_0or1(p,data,bin_vals))
        scores.append( ( normalized_cut(p, inputstates, graph), truthtable ) )
    return sorted(scores, key=lambda x : x[0])


def print_scores(scores):
    for ns,tt in scores:
        print(ns)
        print("\n")
        for t,v in tt.items():
            print("{} : {}".format(t,v))
        print("\n")


def assess_cst_truthtable():
    pass


def test():
    bin_vals = np.array([float(r) for r in range(12)])
    h1 = np.array([0., 1., 3., 5., 2., 0., 0., 0., 0., 0., 0., 0.])
    h2 = np.array([0., 0., 0., 0., 0., 0., 0., 0., 7., 3., 2., 0.])
    h3 = np.array([0., 0., 0., 0., 0., 0., 0., 2., 4., 1., 1., 0.])
    h4 = np.array([0., 0., 0., 0., 1., 1., 2., 8., 3., 0., 0., 0.])
    data = {"00" : h1, "01": h2, "10" : h3, "11" : h4}
    # desired table = 0 1 1 1
    print("\nTest 1: desired table 0 1 1 1, good data\n")
    scores = rank_noncst_tables(data,bin_vals)
    print_scores(scores)

    h1 = np.array([1., 5., 6., 2., 0., 0., 0., 0., 0., 0., 0., 0.])
    h2 = np.array([0., 2., 5., 2., 1., 0., 0., 0., 0., 0., 0., 0.])
    h3 = np.array([0., 0., 0., 1., 3., 7., 2., 1., 0., 0., 0., 0.])
    h4 = np.array([0., 0., 0., 0., 0., 0., 2., 4., 5., 1., 0., 0.])
    data = {"00" : h1, "01": h2, "10" : h3, "11" : h4}
    # desired table = 0 0 1 1
    print("\nTest 2: desired table 0 0 1 1, with state 10 ambiguous\n")
    scores = rank_noncst_tables(data,bin_vals)
    print_scores(scores)


if __name__ == "__main__":
    test()