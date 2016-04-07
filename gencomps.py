# from python standard library
import os
import subprocess
import sys
# third party
import numpy as np
import igraph
# from brainnets
import netgen
import settings
import dataio


def get_global_uw_props(net, props=None):
    """
    For a given network compute the global _unweighted_ properties determined
    by the argument props.

    Parameters
    ----------
    net : igraph.Graph
        the network for which the properties are computed
    props : list
        a list of global unweighted properties as shown in settings.py
        e.g. ["average_clustering", "global_clustering", "assortativity"]
        If props are not given, all global_uw_props listed in settings.py are
        computed.

    Returns
    -------
        measures: a list, where the elements are in the order of the properties
    """
    if props is None:
        props = settings.global_uw_props  # take the default props

    result_dict = {}  # add results in a dict, unwrap later
    if "average_clustering" in props:
        result_dict["average_clustering"] = \
            net.transitivity_avglocal_undirected(mode='zero')
    if "global_clustering" in props:
        result_dict["global_clustering"] = net.transitivity_undirected()
    if "average_path_length" in props:
        result_dict["average_path_length"] = net.average_path_length()
    if "assortativity" in props:
        result_dict["assortativity"] = net.assortativity_degree()
    if "efficiency" in props:
        sp = np.array(net.shortest_paths_dijkstra())
        sp_indices = np.triu_indices_from(sp, 1)
        efficiency_vals = 1. / sp[sp_indices]
        result_dict["efficiency"] = np.average(efficiency_vals)
    if "max_degree" in props:
        result_dict["max_degree"] = max(net.degree(net.vs))
    if "max_kshell" in props:
        result_dict["max_kshell"] = max(net.shell_index())
    return _unwrap(result_dict, props)


def get_global_w_props(net, props=None):
    """
    For a given network compute the global _weighted_ (=correlation) properties
    determined by the argument props.

    Parameters
    ----------
    net : igraph.Graph
        the network for which the properties are computed
    props : list
        a list of global weighted properties as shown in settings.py
        e.g. ["wpl", "wass"]
        If props are not given, all globalWprops listed in settings.py are
        computed.

    Returns
    -------
    measures: a list, where the elements are in the order of the properties
    """
    if props is None:
        props = settings.global_w_props

    result_dict = {}  # add results in a dict, unwrap later
    if "weighted_average_path_length" in props:
        swpls = np.array(
            net.shortest_paths(weights=1. / np.array(net.es["weight"])))
        result_dict["weighted_average_path_length"] = \
            np.average(swpls[np.triu_indices_from(swpls, 1)])
    if "max_strength" in props:
        result_dict["max_strength"] = \
            np.max(net.strength(net.vs, weights=net.es["weight"]))
    if "weighted_clustering" in props:
        result_dict["weighted_clustering"] = \
            net.transitivity_avglocal_undirected(
                mode="zero", weights=net.es["weight"])
    if "weighted_assortativity" in props:
        result_dict["weighted_assortativity"] = net.assortativity(
            net.strength(weights=net.es['weight']), directed=False)
    if 'avg_weight' in props:
        result_dict['avg_weight'] = np.average(np.array(net.es["weight"]))
    return _unwrap(result_dict, props)


def get_global_props_for_density_range(corr_mat, ok_nodes, densities, props,
                                       weighted, include_mst):
    '''
    A wraper to Compute the global (unweighted and weighted) properties
    for a series of network density.

    Parameters
    ----------
    corr_mat : np.array
        2D numpy array with bad nodes.
    ok_nodes : np.array
        the bool blacklist (whitelist)
    density : float
        the network density to use
    props : list
        a list of global properties as shown in settings.py
        If props are not given, all properties listed in settings.py are
        computed.
    weighted : bool
        whether to consider the network as weighted
    include_mst : bool
        whether to include the maximum spanning tree
    '''
    result_dict = {}
    edgeList = netgen.sort_links_by_weight(
        corr_mat, ok_nodes, include_mst=include_mst)
    nNodes = np.sum(ok_nodes)
    nLinksMax = (nNodes * (nNodes - 1)) / 2

    for prop in props:
        result_dict[prop] = np.zeros(len(densities))

    for i, density in enumerate(densities):
        nLinks = int(nLinksMax * density)
        net = netgen.make_net(edgeList[:nLinks], nNodes, weighted=weighted)
        if weighted:
            results = get_global_w_props(net, props)
        else:
            results = get_global_uw_props(net, props)
        for j, prop in enumerate(props):
            result_dict[prop][i] = results[j]
    result_dict[settings.densities_tag] = densities
    return result_dict


def get_node_props(net, props=None):
    """
    Compute the node level properties determined by the argument
    props.

    Parameters
    ----------
    net : igraph.Graph
        the network for which the properties are computed
        the *network should be weighted*
    props : a list of weighted node properties (tags, see :py:mod:settings)
            If props are not given, all node_props listed in
            settings.py are computed.

    Returns
    -------
    measures: a list of lists
        The elements of `measures` correspond to the values of the node
        properties and are provided in the same order as the `props`
        parameter.
        Each list then contains the value of a given property for
        each network node.
    """
    if props is None:
        props = settings.node_props

    result_dict = {}  # add results in a dict, unwrap later

    if "degree" in props:
        result_dict["degree"] = net.degree(net.vs)
    if "strength" in props:
        result_dict["strength"] = net.strength(weights=net.es["weight"])
    if "weighted_betweenness_centrality" in props:
        result_dict["weighted_betweenness_centrality"] = net.betweenness(
            weights=1. / np.array(net.es["weight"]))
    if "betweenness_centrality" in props:
        result_dict["betweenness_centrality"] = net.betweenness()
    if "k_shell" in props:
        result_dict["k_shell"] = net.shell_index()
    if "node_clustering" in props:
        result_dict["node_clustering"] = net.transitivity_local_undirected(
            mode="zero")
    return _unwrap(result_dict, props)


def get_link_props(net, props=None):
    """
    Compute the node level properties determined by the argument
    props.

    Parameters
    ----------
    net : igraph.Graph
        the network for which the properties are computed
        the network should be weighted, if weighted properties
        are computed
    props : an iterable
        a list of node weighted properties as shown in settings.py
        e.g. ["bc"]
        If props are not given, all linkProps listed in
        settings.py are computed.

    Returns
    -------
    resultdict : dict
        a dictionary, where key is the property (e.g. "bc") and
        the value is a np.array containing the value of the
                    measure for each node.
    """
    if props is None:
        props = settings.linkProps
    result_dict = {}  # add results in a dict, unwrap later

    if "edge_betweenness_centrality" in props:
        result_dict["edge_betweenness_centrality"] = net.edge_betweenness()
    if "weighted_edge_betweenness_centrality" in props:
        result_dict["weighted_edge_betweenness_centrality"] = \
            net.edge_betweenness(weights=1. / net.es["weight"])

    return _unwrap(result_dict, props)


def _unwrap(result_dict, props):
    """
    Unwraps the results dict to a list

    Parameters
    ----------
    resultdict : dict
        The dictionary with ``resultdict[prop_key] = value``
    props : list of the properties

    Returns
    -------
    result_list : list
        list of the results
    """
    result_list = []
    for prop in props:
        result_list.append(result_dict[prop])
    return result_list


def get_node_props_from_mat(corr_mat, ok_nodes, density, props, include_mst):
    """
    Get node properties for a specific value of network density

    Parameters
    ----------
    corr_mat : 2D numpy array
        an unfiltered correlation matrix (or equivalent)
    ok_nodes : numpy array, dtype = bool
        numpy bool array, ``ok_nodes[i] = True`` -> node index i is valid
    density : float
        the network density (0.01 = 1%)
    include_mst : bool
        whether or not to include the maximum spanning tree

    Returns
    -------
    result_dict : dict
        dictionary containing the results
        ``result_dict[prop_key] = node_values``
    """
    result_dict = {}
    net = netgen.make_net_from_unfiltered_data(
        corr_mat, ok_nodes, density, include_mst=include_mst, weighted=True)
    results = get_node_props(net, props)
    for i, prop in enumerate(props):
        result_dict[prop] = \
            dataio.expand_1D_node_vals_to_non_blacklisted_array(
                results[i], ok_nodes)
    result_dict[settings.densities_tag] = density
    return result_dict


def get_link_props_from_mat(corr_mat, ok_nodes, density, props, include_mst):
    """
    Get link properties for a specific value of network density
    """
    result_dict = {}
    net = netgen.make_net_from_unfiltered_data(
        corr_mat, ok_nodes, density, include_mst=include_mst, weighted=True)
    results = get_link_props(net, props)
    for i, prop in enumerate(props):
        result_dict[prop] = results[i]
    result_dict[settings.densities_tag] = density
    return result_dict


def comp_link_sim_mat(nets):
    """
    Given a list of networks, computes the link similarity matrix

    Parameters
    ----------
    nets : list
        list of networks

    Returns
    -------
    n_same_links_array : 2D numpy array
        n_same_links_array[i,j] equals to the number of common links
        the two networks have
    """
    n_same_links_array = np.zeros((len(nets), len(nets)))
    for i in range(len(nets)):
        neti = nets[i]
        for j in range(i, len(nets)):
            netj = nets[j]
            netsame = neti.intersection(netj)
            n_same_links_array[i, j] = len(netsame.es)
            n_same_links_array[j, i] = len(netsame.es)
    return n_same_links_array


#
# Module level stuff:
#


def get_best_louvain_partition(graph, weighted, n_it):
    """
    Obtain clusters using the Louvain code by the original authors + Raj's
    randomization.

    (igraph code gives deterministic results for a network)

    Parameters
    ----------
    See :py:func:`get_louvain_partitions`

    Returns
    -------
    resultdict : dict
        dictionary containing the value of modularity and the louvain partition
    """
    retdict = get_louvain_partitions(graph, weighted, n_it)
    modularities = retdict[settings.modularity_tag]
    max_mod_index = np.argmax(modularities)
    return {
        # modularity
        settings.modularity_tag:
        retdict[settings.modularity_tag][max_mod_index],
        # clustering
        settings.louvain_cluster_tag:
        retdict[settings.louvain_cluster_tag][max_mod_index]
    }


def get_louvain_partitions(graph, weighted, n_it):
    """
    Get the n_it number of louvain partitions, for a graph.

    Parameters
    ----------
    graph : igraph.Graph
        the graph for which the modules are computed
    weighted : bool
        should the graph be considered weighted?
    n_it : int
        the number of iterations of the algorithm

    Returns
    -------
    result_dict : dict
        Dictionary containing the level, modularity and clustering
        for each of the iterations.

    Author: Raj, adapted by Rainer
    """
    # put isolated nodes to their own modules:
    components = graph.components()

    # treat isolated nodes separately:
    isolated_nodes = []
    for component in components:
        if len(component) is 1:
            isolated_nodes.append(component[0])

    baseName = "/tmp/" + "louvain_" + str(os.getpid())  # enables parallelism
    # write graph as edge list
    with open(baseName + ".edg", "w") as f:
        if weighted:
            for e in graph.es:
                strToWrite = str(e.source) + " " + \
                    str(e.target) + " " + str(e["weight"]) + "\n"
                f.write(strToWrite)
        else:
            for e in graph.es:
                strToWrite = str(e.source) + " " + str(e.target) + "\n"
                f.write(strToWrite)

    louvainDir = settings.package_dir + "/external_code/gen-louvain/"
    # convert to bin format
    if weighted is True:
        convertRunstr = louvainDir + "convert" + ' -i ' + baseName + ".edg" + \
            ' -o ' + baseName + '.bin' + ' -w ' + baseName + ".weights"
    else:
        convertRunstr = louvainDir + "convert" + ' -i ' + \
            baseName + ".edg" + ' -o ' + baseName + '.bin'
    subprocess.call(convertRunstr.split())

    # Detect the communities

    levels = []
    modularities = []
    nodeLists = []

    for i in range(n_it):
        runstr = louvainDir + "louvain " + baseName + '.bin' + ' -l -1 -q 0'
        if weighted is True:
            runstr += ' -w ' + baseName + ".weights"
        with open(baseName + '.tree', 'w') as fp:
            p = subprocess.Popen(
                runstr.split(), stdout=fp, stderr=subprocess.PIPE)
        _, q = p.communicate()  # modularity score
        modularities.append(float(q))
        # Get the number of hierarchical level
        runstr = louvainDir + "hierarchy " + baseName + ".tree"
        p = subprocess.Popen(
            runstr.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        hLevels, _ = p.communicate()
        hLevel = int(hLevels.split('\n')[0].split(':')[1]) - 1
        levels.append(hLevel)

        # Nodes in the final hierarchical level
        runstr = louvainDir + "hierarchy " + \
            baseName + ".tree" + ' -l ' + str(hLevel)
        p = subprocess.Popen(
            runstr.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        nodeOutput, _ = p.communicate()

        # Read the community
        nodeList = []
        try:
            for line in nodeOutput.split('\n')[:-1]:
                nodeList.append(int(line.split()[1]))
        except:
            print "found..\n"
            print line
            sys.stdout.flush()

        # if isolated nodes are in the end, they are not in the output nodeList

        if len(nodeList) is not len(graph.vs):
            new_max_clu = np.max(nodeList) + 1
            while len(graph.vs) - len(nodeList) is not 0:
                assert len(nodeList) in isolated_nodes, \
                    "something unexpected happened with Louvain"
                nodeList.append(new_max_clu)
                new_max_clu += 1
        nodeLists.append(nodeList)
    return {"levels": levels, settings.modularity_tag: modularities,
            settings.louvain_cluster_tag: nodeLists}


def comp_consensus_partition(clusterings, n_clu_to_be="median"):
    """
    Takes in a number of filtered (i.e. blacklist is removed) partitions
    and computes a consensus cluster using the meta-clustering algorithm.
    (MCLA).

    Parameters
    ----------
    clusterings : list of clusterings / 2D numpy array
    n_clu_to_be : int/str
        How many modules should there (at most) be in the consensus
        partition.
        If "median", the median number of input partitions is used.

    Returns
    -------
    consensus_clu : a numpy array describing the clustering
    """
    import matlab.engine
    eng = matlab.engine.start_matlab()
    eng.addpath(settings.package_dir + "/external_code/ClusterPack-V1.0")
    clusterings = np.array(clusterings) + 1
    clu_nums = []
    for clu in clusterings:
        clu_nums.append(len(np.unique(clu)))
    if n_clu_to_be is "median":
        n_clu_to_be = int(np.median(clu_nums))
    to_matlab = []
    for c in clusterings:
        arr = [int(val) for val in c]
        to_matlab.append(arr)
    matlab_clusterings = matlab.int64(to_matlab)
    consensus_clu = eng.mcla(matlab_clusterings, n_clu_to_be)
    consensus_clu = np.array(consensus_clu[0])-1
    eng.quit()
    return consensus_clu


def comp_partition_sim_mats(
        membership_lists,
        measures=settings.cluster_similarity_measures):
    """
    Computes pairwise clustering similarity measures between different
    partitions (and conditions).

    Parameters
    ----------
    membership_lists : list
        List of membership lists corresponding the partitions
    measures : list
        List of cluster similarity measures

    Returns
    -------
    result_dict : dict
        A dict where the key corresponds to the similarity measure and value
        is a (upper triangular) matrix containing the values of the similarity
        measures.
    """
    membership_lists = np.array(membership_lists, dtype=int)
    n_tot = len(membership_lists)
    result_dict = {}
    for measure in measures:
        sim_mat = np.zeros((n_tot, n_tot))
        for i in range(0, n_tot):
            # print i
            partition1 = membership_lists[i]
            partition1 = [int(partition1[k]) for k in range(len(partition1))]
            for j in range(i, n_tot):
                partition2 = membership_lists[j]
                partition2 = [int(partition2[k])
                              for k in range(len(partition2))]
                sim_mat[i, j] = igraph.compare_communities(
                    partition1, partition2, measure)
                sim_mat[j, i] = sim_mat[i, j]
        result_dict[measure] = sim_mat
    return result_dict


# def computeClusterSimilarityMeasuresNonSymmetric(
#         membership_lists1, membership_lists2,
#         measures=settings.cluster_similarity_measures):
#     """
#     Computes pairwise clustering similarity measures between different
#     clusterings (and conditions).

#     Parameters:
#         List of membership lists corresponding the clusterings

#     Returns:
#         A dict where the key corresponds to the similarity measure and value
#         is a matrix containing the values of the similarity measures.
#     """
#     membership_lists1 = np.array(membership_lists1, dtype=int)
#     membership_lists2 = np.array(membership_lists2, dtype=int)
#     n_tot1 = len(membership_lists1)
#     n_tot2 = len(membership_lists2)
#     result_dict = {}
#     for measure in measures:
#         sim_mat = np.zeros((n_tot1, n_tot2))
#         for i in range(0, n_tot1):
#             clu1 = membership_lists1[i]
#             clu1 = [int(clu1[k]) for k in range(len(clu1))]
#             for j in range(0, n_tot2):
#                 clu2 = membership_lists2[j]
#                 clu2 = [int(clu2[k]) for k in range(len(clu2))]
#                 sim_mat[i, j] = igraph.compare_communities(clu1, clu2,
#                                                           measure)
#         result_dict[measure] = sim_mat
#     return result_dict


def comp_scaled_inclusivity_for_ref_partition(
        ref_partition, other_partitions, normalize=False):
    """
    Computes the SI-measure for all nodes wrt. to the ref_partition

    See `Assessing the consistency of community structure in complex networks
    <http://pre.aps.org/abstract/PRE/v84/i1/e016111>`_ for more information.

    For a node i, the SI value equals to

    .. math::

        SI(i) = \sum_{n} \\frac{|R_i \cap C_n^i|^2}{|R^i| |C_n^i|}

    Where the summation is over the set of other partitions :math:`\\{C_n\\}`.
    (:math:`|C_n^i|` denotes the number of nodes in the module `C_n^i` )

    Parameters
    ----------
    ref_partition : a np.array
        numpy array with shape (n_nodes,) containing the reference membership
        list
    other_partitions : list of numpy arrays / 2D numpy array
        iterable containing the other partitions as membership lists
    normalize : bool
        whether to normalize by the number of comparisons made, so that
        nodewise SI values are in range [0,1]

    Returns
    -------
    node_SIs : numpy array
        the node-wise SI values

    See also
    --------
    comp_scaled_inclusivity : SI for a set of partitions
    """
    assert isinstance(
        ref_partition, np.ndarray), "ref_partition is not a numpy array"
    assert isinstance(other_partitions[0], np.ndarray), \
        "argument other_partitions should be an iterable of numpy arrays"
    node_SIs = np.zeros(len(ref_partition))
    for i, label in enumerate(ref_partition):
        print i
        ref_community = (ref_partition == label)  # np.array of bool elements
        size_ref_community = np.sum(ref_community)
        for j in range(len(other_partitions)):
            label_other_community = other_partitions[j][i]
            other_community = (other_partitions[j] == label_other_community)
            intersection_size = np.sum(ref_community * other_community)
            size_other_community = np.sum(other_community)
            node_SIs[i] += (intersection_size * intersection_size) / \
                float(size_other_community * size_ref_community)
            assert size_ref_community >= 1
    if normalize:
        return node_SIs / float(len(other_partitions))
    else:
        return node_SIs


def comp_scaled_inclusivity(partitions, normalize=True):
    """
    Computes the SI-measure for all nodes between partitions.

    For a node i, the SI value equals to

    .. math::

        SI(i) = \sum_{n, m} \\frac{|C_n^i \cap C_m^i|^2}{|C_n^i| |C_m^i|}

    Where the summation is over the set of all partition pairs :math:`C_n,C_m`.
    (:math:`|C_n^i|` denotes the number of nodes in the module `C_n^i` )

    See `Assessing the consistency of community structure in complex networks
    <http://pre.aps.org/abstract/PRE/v84/i1/e016111>`_ for more information.

    Parameters
    ----------
    partitions : (list of numpy arrays or a 2D np.ndarray)
        list of membership lists describing the different partitions
    normalize : bool
        whether to normalize by the number of comparisons made, so that
        nodewise SI values lie within range [0,1]

    Returns
    -------
    node_SIs : np.array
        the nodewise SI values

    See Also
    --------
    comp_scaled_inclusivity_for_ref_partition : SI with respect to a certain
        partition
    """
    assert isinstance(
        partitions[0], np.ndarray), "partitions are not numpy arrays"
    # no nans allowed:
    assert np.logical_not(np.isnan(partitions)).all(), \
        "the input partitions should not contain nan values"
    n_nodes = len(partitions[0])
    n_partitions = len(partitions)
    node_SIs = np.zeros(n_nodes)
    for i in range(0, n_partitions):
        print "partition ", i, "of in total ", n_partitions
        partition_i = partitions[i]
        uniques_i = np.unique(partition_i)
        for j in range(i + 1, n_partitions):
            partition_j = partitions[j]
            uniques_j = np.unique(partition_j)
            partition_sim_mat = {}
            for ilabel in np.sort(uniques_i):
                ilabelclu = (partition_i == ilabel)
                ilabelclusize = np.sum(ilabelclu)
                for jlabel in np.sort(uniques_j):
                    jlabelclu = (partition_j == jlabel)
                    jlabelclusize = np.sum(jlabelclu)
                    intersection_size = np.sum(ilabelclu * jlabelclu)
                    partition_sim_mat[
                        (ilabel, jlabel)] = (intersection_size ** 2 /
                                             (float(ilabelclusize *
                                                    jlabelclusize)))
                    # print (ilabel, jlabel), \
                    #      partition_sim_mat[(ilabel, jlabel)]
            for node in range(n_nodes):
                clu_i = (partition_i[node])
                clu_j = (partition_j[node])
                node_SIs[node] += partition_sim_mat[clu_i, clu_j]
    if normalize:
        node_SIs = node_SIs / (float(n_partitions * (n_partitions - 1) * 0.5))
    return node_SIs

# def matchClustersHungarianAlgo(cluster1, cluster2):
#     """
#     Match the clusters using Hungarian Algorithm
#     """
#     from rpy2.robjects import r
#     from rpy2.robjects import IntVector
#     import rpy2.robjects.numpy2ri
#     rpy2.robjects.numpy2ri.activate()
#     # cluster labels should be from 1:n
#     newCluster1 = np.zeros(len(cluster1))
#     newCluster2 = np.zeros(len(cluster2))
#     for i, label in enumerate(np.unique(cluster1)):
#         newCluster1[cluster1 == label] = i + 1
#     for i, label in enumerate(np.unique(cluster2)):
#         newCluster2[cluster2 == label] = i + 1
#
#     r.source(settings.package_dir + 'rfiles/community_structure.R')
#     hungarianmatch = r["hungarianmatch"]
#     cluster = hungarianmatch(IntVector(newCluster1), IntVector(newCluster2))
#     cluster = map(int, cluster)
#     result = np.array([c for c in cluster])
#     return result
