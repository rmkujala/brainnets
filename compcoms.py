# import gc
# third party
import numpy as np
import igraph
# brainnets imports
import aux
import config
import comp_helpers as ch
import dataio
import fname_conventions as fnc
import gencomps
import netgen
import settings


def tag_to_igraph_comdet_method(tag):
    tag_to_method = {f: m for f, m in vars(igraph.Graph).items()
                     if (f == "community_" + tag)}
    try:
        return tag_to_method["community_" + tag]
    except KeyError as e:
        raise e


def igraph_com_det_method_to_tag(igraph_com_det_method):
    assert igraph_com_det_method.__name__[:10] == "community_"
    return igraph_com_det_method.__name__[10:]


def comp_communities_igraph(cfg, com_det_method, com_det_options_dict=None):
    """
    Computes communities for a certain network density for all
    correlation matrices in ``cfg['all_fnames']``
    The results are saved to the output folder (``cfg['outdata_dir']``)

    Parameters
    ----------
    cfg : a brainnets config dict
    com_det_method: str, or igraph function returning
    com_det_options_dict: dict

    Returns
    -------
    coms :
        the communities as a membership list
    """
    if isinstance(com_det_method, str):
        com_det_method = tag_to_igraph_comdet_method(com_det_method)
    if com_det_options_dict is None:
        com_det_options_dict = {}
    config.require(cfg, ['all_fnames',
                         'blacklist_fname',
                         'density',
                         'include_mst',
                         'n_cpus',
                         'n_it_comdet'
                         ]
                    )
    arg_list = [(fname, cfg, com_det_method, com_det_options_dict)
               for fname in cfg['all_fnames']]
    coms = ch.run_in_parallel(_compute_coms_worker, arg_list, cfg['n_cpus'])
    return coms


def _compute_coms_worker(args):
    """
    Computes Louvain communities.
    """
    fname, cfg, com_det_method, com_det_options= args
    coms = []
    graph = netgen.get_graph_from_bare_data(
        fname, cfg['blacklist_fname'], cfg['density'],
        include_mst=cfg['include_mst'], weighted=False)
    membershiplists = []
    for i in range(cfg['n_it_comdet']):
        clustering = com_det_method(graph, **com_det_options)
        if isinstance(clustering, igraph.clustering.VertexDendrogram):
            clustering = clustering.as_clustering()
        membershiplists.append(clustering.membership)
    coms = np.array(membershiplists)
    ok_nodes = dataio.get_ok_nodes(cfg['blacklist_fname'])
    # expand communities to non-filtered indices
    unfiltered_coms = []
    for i, com in enumerate(coms):
        uf_com = dataio.expand_1D_node_vals_to_non_blacklisted_array(
            com, ok_nodes
        )
        unfiltered_coms.append(uf_com)
    unfiltered_coms = np.array(unfiltered_coms)
    com_det_method_tag = igraph_com_det_method_to_tag(com_det_method)
    out_fname = fnc.get_ind_fname(fname, cfg, com_det_method_tag)
    out_dict = {com_det_method_tag : unfiltered_coms,
                settings.config_tag         : cfg}
    dataio.save_pickle(out_fname, out_dict)
    print "finished " + fname
    return unfiltered_coms


def comp_louvain_communities(cfg):
    """
    Computes louvain communities for a certain network density for all
    correlation matrices in ``cfg['all_fnames']``
    The results are saved to the output folder (``cfg['outdata_dir']``)

    Currently only the unweighted louvain method is used.

    Parameters
    ----------
    cfg : a brainnets config dict

    Returns
    -------
    coms :
        the communities as a membershiplist
    mods : numpy array
        the corresponding values of modularity
    """
    config.require(cfg, ['all_fnames', 'blacklist_fname', 'density',
                         'n_it_comdet', 'include_mst', 'n_cpus'])
    argList = [(fname, cfg) for fname in cfg['all_fnames']]
    comsAndModularities = ch.run_in_parallel(
        _compute_louvain_coms_worker, argList, cfg['n_cpus'])
    coms = [comsAndModularities[i][0] for i in range(len(comsAndModularities))]
    mods = [comsAndModularities[i][1] for i in range(len(comsAndModularities))]
    return coms, mods


def _compute_louvain_coms_worker(args):
    """
    Computes Louvain communities.
    """
    fname, cfg = args
    coms = []
    mods = []
    print "started " + fname
    graph = netgen.get_graph_from_bare_data(
        fname, cfg['blacklist_fname'], cfg['density'],
        include_mst=cfg['include_mst'], weighted=False)
    louvain_coms_dict = \
        gencomps.get_louvain_partitions(graph, False, cfg['n_it_comdet'])
    coms.extend(louvain_coms_dict[settings.louvain_cluster_tag])
    coms = np.array(coms)
    ok_nodes = dataio.get_ok_nodes(cfg['blacklist_fname'])
    # expand communities to non-filtered indices
    unfiltered_coms = []
    for i, com in enumerate(coms):
        uf_com = dataio.expand_1D_node_vals_to_non_blacklisted_array(
            com, ok_nodes
        )
        unfiltered_coms.append(uf_com)
    unfiltered_coms = np.array(unfiltered_coms)
    mods.extend(louvain_coms_dict[settings.modularity_tag])
    mods = np.array(mods)
    out_fname = fnc.get_ind_fname(fname, cfg, settings.louvain_cluster_tag)
    out_dict = {settings.louvain_cluster_tag: unfiltered_coms,
                settings.modularity_tag: mods,
                settings.config_tag: cfg}
    dataio.save_pickle(out_fname, out_dict)

    print "finished " + fname
    return unfiltered_coms, mods


def comp_consensus_partition(cfg, fnames_tag, out_fname,
                             n_clu_for_mcla='median',
                             n_to_consider=None,
                             comdet_tag=None):
    """
    Computes a consensus partition.

    Parameters
    ----------
    cfg : dict
        a brainnets config dictionary
    fnames_tag : str
        the filenames group for which the consensus partition is
        computed
    out_fname : str
        the filename to which the consensus partition is stored
    n_clu_for_mcla : int or "median"
        maximum number or clusters in the consensus partition
        if "median", the median number is used as the max number
        of clusters in the consensus partition
    n_to_consider : int/str, optional
        number of partitions to consider for obtaining consensus
        defaults to considering _all_ partitions
        if "best" uses the partition with maximum modularity
        if available
    comdet_tag: str, optional
        e.g. "infomap"
        defaulting to settings.louvain_cluster_tag (legacy)

    Returns
    -------
    out_dict : dict
        dictionary containing the consensus partition
    """
    config.require(cfg, [fnames_tag, 'blacklist_fname', 'density'])

    ok_nodes = dataio.get_ok_nodes(cfg['blacklist_fname'])
    if comdet_tag is None:
        comdet_tag = settings.louvain_cluster_tag

    # load clusterings
    clusterings = None
    ok_nodes = dataio.get_ok_nodes(cfg['blacklist_fname'])
    for fname in cfg[fnames_tag]:
        indfname = fnc.get_ind_fname(fname, cfg, comdet_tag)
        data = dataio.load_pickle(indfname)
        clus_raw = data[comdet_tag]

        assert len(clus_raw[0]) >= np.sum(ok_nodes)
        if n_to_consider is not None:
            if isinstance(n_to_consider, int):
                clus_raw = clus_raw[:n_to_consider]
            elif n_to_consider == 'best':
                max_mod_i = np.argmax(data[settings.modularity_tag])
                clus_raw = clus_raw[max_mod_i]
                clus_raw = clus_raw.reshape(1, len(clus_raw))
            else:
                assert isinstance(n_to_consider, int) or n_to_consider == 'best', \
                    "n_to_consider should be an integer!"

        clus = clus_raw[:, ok_nodes]
        if clusterings is None:
            # for first encounter
            clusterings = np.copy(clus)
        else:
            clusterings = np.vstack((clusterings, clus))

    # this should hold usually, unless you have a non-standard workflow:
    # (added for making sure a bug does not exist anymore)
    assert len(clusterings) == len(clus) * len(cfg[fnames_tag])

    # print len(clusterings), n_clu_for_mcla
    consensus_clu = gencomps.comp_consensus_partition(
        clusterings, n_clu_for_mcla)
    consensus_clu = dataio.expand_1D_node_vals_to_non_blacklisted_array(
        consensus_clu, ok_nodes, default_value=-1)
    out_dict = {comdet_tag: consensus_clu,
                settings.config_tag: cfg}

    dataio.save_pickle(out_fname, out_dict)
    return out_dict


def comp_scaled_inclusivity_for_two_fname_groups(cfg):
    config.require(
        cfg, ["density", "group_1_mat_fnames", "group_2_mat_fnames"])
    fname_groups = [cfg['group_1_mat_fnames'], cfg['group_2_mat_fnames']]
    for i, fname_group in enumerate(fname_groups):
        clus = []
        for mat_fname in fname_group:
            clusters_fname = fnc.get_ind_fname(
                mat_fname,
                cfg,
                settings.louvain_cluster_tag
            )
            subject_clusters = dataio.load_pickle(clusters_fname)
            clus.append(subject_clusters[settings.louvain_cluster_tag])
        partitions = aux.expand_first_axis(np.array(clus))
        partitions = partitions[:, dataio.get_ok_nodes(cfg['blacklist_fname'])]
        assert np.logical_not(np.isnan(partitions)).all()
        node_SIs = gencomps.comp_scaled_inclusivity(partitions)
        out_dict = {settings.scaled_inclusivity_tag:
                    node_SIs, settings.config_tag: cfg}
        out_fname = fnc.get_group_fname(
            cfg, settings.scaled_inclusivity_tag, i)
        dataio.save_pickle(out_fname, out_dict)


def comp_consensus_scaled_inclusivity(cfg, group_id, n_to_consider=None):
    """
    Parameters
    ----------
    cfg : dict
        brainnets config dictionary
    group_id : int
        0 or 1 -- the group for which the scaled inclusivity should be computed
    """
    config.require(
        cfg, ["density", "group_1_mat_fnames", "group_2_mat_fnames"])

    if group_id == 0:
        fname_group = cfg['group_1_mat_fnames']

    elif group_id == 1:
        fname_group = cfg['group_2_mat_fnames']
    else:
        raise Error('Param group_id should be either 0 or 1')
    consenus_com_fname = fnc.get_group_fname(
        cfg, settings.louvain_consensus_tag, group_id)
    consensus_com = \
        dataio.load_pickle(consenus_com_fname)[settings.louvain_cluster_tag]

    clus = []
    for mat_fname in fname_group:
        clusters_fname = fnc.get_ind_fname(
            mat_fname,
            cfg,
            settings.louvain_cluster_tag
        )
        data = dataio.load_pickle(clusters_fname)
        subject_clusters = data[settings.louvain_cluster_tag]

        if n_to_consider is not None:
            if isinstance(n_to_consider, int):
                subject_clusters = subject_clusters[:n_to_consider]
            elif n_to_consider == 'best':
                max_mod_i = np.argmax(data[settings.modularity_tag])
                subject_clusters = subject_clusters[max_mod_i]
                subject_clusters = subject_clusters.reshape(
                    1, len(subject_clusters))
            else:
                assert isinstance(n_to_consider, int) or n_to_consider == 'best', \
                    "n_to_consider should be an integer!"
        clus.append(subject_clusters)

    partitions = aux.expand_first_axis(np.array(clus))
    ok_nodes = dataio.get_ok_nodes(cfg['blacklist_fname'])
    partitions = partitions[:, ok_nodes]
    consensus_com = consensus_com[ok_nodes]
    assert np.logical_not(np.isnan(partitions)).all()
    assert len(consensus_com) == len(partitions[0])

    node_SIs = gencomps.comp_scaled_inclusivity_for_ref_partition(
        consensus_com, partitions, normalize=True)
    out_dict = {settings.scaled_inclusivity_tag:
                node_SIs, settings.config_tag: cfg}
    out_fname = fnc.get_group_fname(
        cfg, settings.louvain_consensus_si_tag, group_id)
    dataio.save_pickle(out_fname, out_dict)
