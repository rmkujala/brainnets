from brainnets import settings
import numpy as np
import igraph

# def comp_stable_community_adj_mat_based_on_best(fnames, density,
#                                                 blacklist_fname):
#     okNodes = dataio.get_ok_nodes(blacklist_fname)
#     communities = []
#     for fName in fnames:
#         com = dataio.loadPickle(
#                dataio.getLouvainClusteringIndividualFileName(
#                       fName, density))[settings.louvain_cluster_tag][okNodes]
#         communities.append(com)
#     adjMat = make_adj_mat_from_partitions(communities, normalize=True)
#     return adjMat


def make_adj_mat_from_partitions(partitions, normalize=True):
    """
    Given a number of network partitions computes the weight
    matrix, where the weight of a link reflects the number of times
    the two nodes are in the same community.


    Parameters
    ----------
    partitions: list of lists or numpy array
        the network partitions (a list of membership lists)
    normalize: bool
        whether to normalize the weights to the range [0,1]
    """
    partitions = np.array(partitions)
    n_nodes = len(partitions[0])
    adj_mat = np.zeros((n_nodes, n_nodes))
    for i, com in enumerate(partitions):
        # print i, for logging
        print i
        com_labels = np.unique(com)
        for com_label in com_labels:
            # indices of common nodes
            com_nodes = np.nonzero(com == com_label)[0]
            # cartesian product.. hope this works!
            com_size = len(com_nodes)
            coords = (
                [np.tile(com_nodes, com_size), np.repeat(com_nodes, com_size)])
            adj_mat[coords] += 1

#            for i in range(len(comNodes)):
#                comNodeI = comNodes[i]
#                for j in range(i+1, len(comNodes)):
#                    comNodeJ = comNodes[j]
#                    adj_mat[comNodeI, comNodeJ] += 1
    # adj_mat += adj_mat.T #make it symmetric
    if normalize:
        adj_mat = adj_mat / float(len(partitions))
    return adj_mat


def get_module_sizes(clu):
    """
    Gives the module sizes as a numpy vector.

    Parameters
    ----------
    clu : numpy vector
        membership list (blacklist filtered)

    Returns
    -------
    sizes : numpy array
        the sizes of the modules
    """
    sorted_clu_labels = np.sort(np.unique(clu))
    sizes = np.zeros(len(sorted_clu_labels))
    for i, cluLabel in enumerate(sorted_clu_labels):
        sizes[i] = np.sum((clu == cluLabel))
    return sizes


def get_module_sizes_as_dict(clus):
    """
    Gives the module sizes as a dict

    Parameters
    ----------
    clu : numpy vector
        membership list (blacklist filtered)

    Returns
    -------
    size_dict : dict
        the sizes of the modules
        ``size_dict[module_index] = module_size``
    """
    size_dict = {}
    for i in np.unique(clus):
        size_dict[i] = np.sum(clus == i)
    return size_dict


def print_cluster_sizes(partition):
    """
    Prints cluster sizes to standard output.
    The module with the cluster label:
        :py:attr:`settings.undef_clu_label`
    Will be ignored.

    Parameters
    ----------
    partition : list or numpy array
        The partition, which contains the modules.
        (A membership list)
    """
    clu_labels = np.setdiff1d(np.unique(partition), [settings.undef_clu_label])
    for label in clu_labels:
        print np.sum(partition == label),
    print ""


def get_stable_communities(adj_mat, alpha):
    """
    Returns the core communities as a cluster assingment list
    (
    stable community cores:
    http://link.springer.com/chapter/10.1007%2F978-3-642-30287-9_10#page-1
    )
    (starting from zero - the simplest way to present a clustering)
    e.g.::

        [0, 0, 1, 1, 2, 2, 0, 1, 2]

    Parameters
    ----------
    adj_mat :
        numpy.ndarray, the symmetric adjacency matrix
        (only the lower diag used though)
    alpha :
        float
        the treshold for the stable "graph" / matrix

    Returns
    -------
    components_as_array : numpy array
        As above ::

            [0, 0, 1, 1, 2, 2, 0, 1, 2])

    components : list
        The commponents as a list of lists::

            [[0, 1, 6], [2, 3, 7], [4, 5, 8]]

    See also
    --------
    mask_small_communities : to get rid of the small non-stable communities
    """
    n = len(adj_mat)
    print np.sum(adj_mat > alpha) / (n * (n - 1.))
    graph = igraph.Graph(n)
    for i in range(n):
        for j in range(i, n):
            if adj_mat[i, j] > alpha:
                graph.add_edge(i, j)
    components = graph.components()
    components_as_array = np.zeros(n, dtype=int)
    components_as_list = []
    for i, com in enumerate(components):
        components_as_array[com] = i
        components_as_list.append(com)
    return components_as_array, components_as_list


def mask_small_communities(partition, clu_size_threshold, undef_clu_label=-1):
    """
    Removes clusters of size < clu_size_threshold to one big cluster marked
    by the param undef_clu_label

    Parameters
    ----------
    partition : list or numpy array
        a cluster membership list
    clu_size_threshold : int or float
        clusters smaller than ``clu_size_threshold`` that will be grouped
        to one big one

    Returns
    -------
    partition :
        the filtered partition
    """
    partition = np.array(partition)
    clu_labels = np.unique(partition)
    for cluLabel in clu_labels:
        clu = partition == cluLabel
        if np.sum(clu) < clu_size_threshold:
            partition[clu] = undef_clu_label
    return partition


def comp_cluster_adj_mat(partition, net):
    """
    Computes the clusterwise adj - matrix(with 'self loops').

    Parameters
    ----------
    partition : list or numpy array
        a membership list:
            partition[nodeindex] = clulabel
            partitions should be indexed from zero,
            Module assignments with :py:attr:`settings.undef_clu_label`
            are neglected.
    net : igraph.Graph
        The original links are taken from here.
        (e.g. ``net.es[0].source`` equals to the "link 0 source node")

    Returns
    -------
    clu_adj_mat :
        numpy array
        The computed cluster adjacency matrix
    """
    n_clus = len(np.setdiff1d(partition, [settings.undef_clu_label]))
    clu_adj_mat = np.zeros((n_clus, n_clus))
    for e in net.es:
        # have to loop over each as += operation does not work
        source_clu = partition[e.source]
        target_clu = partition[e.target]
        assert (source_clu >= 0 and target_clu >= 0)
        clu_adj_mat[source_clu, target_clu] += 1
        if source_clu != target_clu:
            clu_adj_mat[target_clu, source_clu] += 1
    return clu_adj_mat



def make_zero_indexed_clustering(clustering,
                                 undef_clu_label=settings.undef_clu_label):
    """
    Takes a clustering ::

        [1, 1, 1, 1, 3, 2, 3, 1, ..., 11, 100])

    and relabels the clusters starting from zero ::

        [0, 0, 0, 0, 2, 1, 2, 0, ..., 10, 99]

    Parameters
    ----------
    clustering : list or numpy array
        membership list
    undef_clu_label : int
        The label for nodes, for which clustering is not defined.
        Defaults to:
            :py:attr:`settings.undef_clu_label`

    Returns
    -------
    clustering_new : numpy array
        The modified clustering
    """
    real_clu_labels = np.sort(
        np.setdiff1d(np.unique(clustering), [undef_clu_label]))
    clustering_new = np.zeros(np.shape(clustering), dtype=np.int64)
    if undef_clu_label is not None:
        bad_indices = ((np.array(clustering) == undef_clu_label) + np.isnan(clustering))
        clustering_new[bad_indices] = undef_clu_label
    for i, cluLabel in enumerate(real_clu_labels):
        clustering_new[clustering == cluLabel] = i
    return clustering_new


def match_clusters_greedy(clus1, clus2):
    """
    Match partitions clus1, and clus2 reasonably using a greedy clustering
    scheme. Matches the partition with more clusters to the one with less
    clusters. If equal number of clusters, clus2 will be matched to clus1
    The matching could probably be upgraded based on greedy ribbon
    overlap maximisation.

    Parameters
    ----------
        clus1, clus2 : list
            cluster assignments of clu1 [0,1,9,8,6,10,9,8,...]
    Returns
    -------
        new_clus_1, newClus2 : the new mathed clusters
    """
    clus1 = np.array(clus1)
    clus2 = np.array(clus2)
    clusizes_dict_1 = get_module_sizes_as_dict(clus1)
    clusizes_dict_2 = get_module_sizes_as_dict(clus2)
    n_clus_1 = len(clusizes_dict_1)
    n_clus_2 = len(clusizes_dict_2)
    if n_clus_1 <= n_clus_2:
        refClus = clus1
        refCluSizesDict = clusizes_dict_1
        matchClus = clus2
        matchCluSizesDict = clusizes_dict_2
    else:
        refClus = clus2
        refCluSizesDict = clusizes_dict_2
        matchClus = clus1
        matchCluSizesDict = clusizes_dict_1

#    for matchCluLabel, _ in sorted(matchCluSizesDict.iteritems(),
    # key=lambda x: -x[1]):
#        matchClu = (matchClus == matchCluLabel)
#        overLapDict = {}
#        for refCluLabel in np.unique(refClus):
#            refClu = (refClus==refCluLabel)
#            overLap = refCluFreeDict[refCluLabel]*np.sum(refClu*matchClu)
#            overLapDict[refCluLabel] = overLap
# get largest overlap
#        max_key, max_val = sorted(overLapDict.iteritems(),
        # key=lambda x: -x[1])[0]
#        if max_val > 0:
#            matchClusNew[matchClu] = max_key
#            refCluFreeDict[max_key] = 0
#        else:
# find
#            freeKey = np.min(otherFreeIndices)
#            otherFreeIndices[otherFreeIndices==freeKey] = \
        # np.max(otherFreeIndices)+1
#            matchClusNew[matchClu] = freeKey
# pass # do nothing as no reasonable match was found
    # compute overlap matrix
    overLapMat = np.zeros((len(refCluSizesDict), len(matchCluSizesDict)))
    origRefCluLabels = np.array(refCluSizesDict.keys())
    origMatchCluLabels = np.array(matchCluSizesDict.keys())
    for i, refCluLabel in enumerate(refCluSizesDict.keys()):
        refClu = (refClus == refCluLabel)
        for j, matchCluLabel in enumerate(matchCluSizesDict.keys()):
            matchClu = (matchClus == matchCluLabel)
            overLapMat[i, j] = np.sum(refClu * matchClu)
    # look up for free indices"
    refCluFreeDict = {}
    otherFreeIndices = []
    for key in range(np.min(refClus), np.max(refClus) + np.max(matchClus) +
                     np.min(matchClus)):
        if key in refClus:
            refCluFreeDict[key] = 1
        else:
            otherFreeIndices.append(key)  # not used but free
    otherFreeIndices = np.array(otherFreeIndices)

    matchClusNew = np.copy(matchClus)
    # search for largest i, j in over overLapMat

    origMatchCluLabelsRemaining = set(origMatchCluLabels)
    while np.sum(overLapMat) != 0:
        i, j = np.unravel_index(overLapMat.argmax(), overLapMat.shape)
        refCluLabel = origRefCluLabels[i]
        matchCluLabel = origMatchCluLabels[j]
        matchClu = (matchClus == matchCluLabel)
        matchClusNew[matchClu] = refCluLabel
        overLapMat[i, :] = 0
        overLapMat[:, j] = 0
        refCluFreeDict[refCluLabel] = 0
        origMatchCluLabelsRemaining.remove(matchCluLabel)
    while len(origMatchCluLabelsRemaining) > 0:
        matchCluLabel = origMatchCluLabelsRemaining.pop()
        matchClu = (matchClus == matchCluLabel)
        # find
        freeKey = np.min(otherFreeIndices)
        otherFreeIndices[otherFreeIndices == freeKey] = np.max(
            otherFreeIndices) + 1
        matchClusNew[matchClu] = freeKey

    if n_clus_1 <= n_clus_2:
        assert n_clus_1 == len(np.unique(refClus))
        assert n_clus_2 == len(np.unique(matchClusNew))
        return refClus, matchClusNew
    else:
        assert n_clus_1 == len(np.unique(matchClusNew))
        assert n_clus_2 == len(np.unique(refClus))
        return matchClusNew, refClus
