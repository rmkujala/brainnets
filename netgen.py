# -*- coding: utf-8 -*-
"""
This module contains functions for thresholding matrices
and outputting links/networks.
"""
import numpy as np
import igraph
import dataio


def get_graph_from_bare_data(corr_mat_fname, blacklist_fname, density,
                             include_mst=False, weighted=False):
    """
    Extracts a graph from raw data.

    Parameters
    ----------
    corr_mat_fname : str
        path to the file containing the correlation matrix.
    blacklist_fname : str
        path to the bool blacklist
    density : float
        the network density to use
    include_mst : bool
        whether to include the maximum spanning tree
    weighted : bool
        whether to consider the network as weighted

    Returns
    -------
    net : igraph.Graph
        the network
    """
    corr_mat = dataio.load_adj_matrix_from_mat(corr_mat_fname)
    ok_nodes = dataio.get_ok_nodes(blacklist_fname)
    net = make_net_from_unfiltered_data(
        corr_mat,
        ok_nodes,
        density,
        include_mst=include_mst,
        weighted=weighted)
    return net


def _get_filtered_triu_adj_mat_copy(matrix, ok_nodes):
    """
    Takes only the nodes listed in ok_nodes into account.

    Parameters
    ----------
    matrix : np.array
        2D matrix with bad nodes
    ok_nodes : numpy bool array

    Returns
    -------
    m : np.array
        a copy of the matrix where the bad nodes have been removed
    """
    m = matrix.copy()
    m = m[ok_nodes, :]
    m = m[:, ok_nodes]
    return np.triu(m, 1)


def make_net_from_unfiltered_data(corr_mat, ok_nodes, density, include_mst=False,
                                  weighted=False):
    """
    Constructs a net from unfiltered data.

    Parameters
    ----------
    corr_mat : np.array
        2D numpy array with bad nodes.
    ok_nodes : np.array
        the bool blacklist (whitelist)
    density : float
        the network density to use
    include_mst : bool
        whether to include the maximum spanning tree
    weighted : bool
        whether to consider the network as weighted

    Returns
    -------
    net : igraph.Graph
    """
    assert 0 <= density <= 1
    edgelist = sort_links_by_weight(corr_mat, ok_nodes, include_mst)

    nNodes = sum(ok_nodes)
    nLinksMax = (nNodes * (nNodes - 1)) / 2
    nLinks = int(nLinksMax * density)
    edgelist = edgelist[:nLinks]

    return make_net(edgelist, nNodes, weighted)


def get_treshold_value(corr_mat, ok_nodes, density, include_mst=False):
    """
    Constructs a net from unfiltered data.

    Parameters
    ----------
    corr_mat : np.array
        2D numpy array with bad nodes.
    ok_nodes : np.array
        the bool blacklist (whitelist)
    density : float
        the network density to use
    include_mst : bool
        whether to include the maximum spanning tree

    Returns
    -------
    threshold: float
        the weight corresponding to the last considered link
        (i.e. no threshold)
    """
    assert 0 <= density <= 1
    edgelist = sort_links_by_weight(corr_mat, ok_nodes, include_mst)
    n_nodes = sum(ok_nodes)
    n_links_max = (n_nodes * (n_nodes - 1)) / 2
    n_links = int(n_links_max * density)
    return edgelist[n_links]['weight']


def make_net(edgelist, nNodes, weighted):
    '''
    Create the network given the edgelist and number of nodes

    Parameters
    ----------
    weighted : (boolean)
        Whether weights are to be considered or not
    '''
    graph = igraph.Graph(nNodes)

    graph.add_edges(zip(edgelist['node1'], edgelist['node2']))
    if weighted is True:
        # graph.es['weight'] = 1
        graph.es['weight'] = edgelist['weight']

    # for n1, n2, w in edgelist:
    #     graph[n1, n2] = w
    return graph


def make_full_weighted_net_from_weight_mat(matrix, ok_nodes, return_weights=False):
    """
    Takes in an adjacency/correlation matrix, and constructs an undirected
    weighted network
    """
    nNodes = np.sum(ok_nodes)
    graph = igraph.Graph(nNodes)

    triu_indices = np.triu_indices_from(matrix, 1)
    edgelist = np.array(triu_indices).T
    graph.add_edges(edgelist)
    weights = matrix[triu_indices]
    graph.es["weight"] = weights
    if return_weights:
        return graph, weights
    return graph


def sort_links_by_weight(corr_mat, ok_nodes, include_mst):
    """
    Sort the links by their link-weight

    Parameters
    ----------
    corr_mat : np.array
        2D numpy array with bad nodes.
    ok_nodes : np.array
        the bool blacklist (whitelist)
    include_mst : Bool
        If true add the maximum spanning tree to the begining of sorted list

    Returns
    -------
    edgelist : numpy structrued array (node1, node2, weight)
        array([(0, 1, 1.0), (0, 3, 0.5), (2, 3, 0.5), (0, 4, 0.7), (1, 4, 0.4)],
              dtype=[('node1', '<i4'), ('node2', '<i4'), ('weight', '<f8')])
    """
    up_diag_matrix = _get_filtered_triu_adj_mat_copy(corr_mat, ok_nodes)
    n = len(up_diag_matrix)
    minVal = np.min(up_diag_matrix)
    minValMinusOne = np.min(up_diag_matrix) - 1
    # So that possible overflows don't go unnoticed
    assert minValMinusOne < minVal

    initEdges = np.array(np.triu_indices_from(up_diag_matrix, 1)).T
    weights = up_diag_matrix[np.triu_indices_from(up_diag_matrix, 1)]
    nLinksMax = (n * (n - 1)) / 2
    nLinksMST = 0
    edgelist = np.zeros(
        nLinksMax, dtype=[('node1', 'i4'), ('node2', 'i4'), ('weight', 'f8')])

    # Get the maximum spanning tree (Multyplying the weights by -1 does the
    # trick)
    if include_mst:
        g = igraph.Graph(n, list(initEdges), directed=False)
        mst = g.spanning_tree(-1 * weights, return_tree=False)
        for i, ei in enumerate(mst):
            edge = g.es[ei]
            edgelist[i] = edge.source, edge.target, weights[ei]
            # Take these links away from the orig. mat
            up_diag_matrix[edge.source, edge.target] = minValMinusOne
        nLinksMST = len(mst)

    # How many links we still need to take after (possible) MST:
    nLinksYetToTake = np.max([nLinksMax - nLinksMST, 0])  # mst already there

    # Get the next largest indices
    up_diag_matrix[np.tril_indices_from(up_diag_matrix, 0)] = minValMinusOne
    mflat = up_diag_matrix.flatten()
    flatindices = mflat.argsort()[::-1][:nLinksYetToTake]
    edgelist[nLinksMST:]['node1'], edgelist[nLinksMST:][
        'node2'] = np.unravel_index(flatindices, (n, n))
    edgelist[nLinksMST:]['weight'] = mflat[flatindices]

    return edgelist
