# brainnets:
import numpy as np
import config
import comp_helpers as ch
import dataio
import netgen
import settings
import fname_conventions as fnc


def comp_link_distances(cfg):
    """
    Computes all link distances (in MNI space) over a range of network
    densities.

    Parameters
    ----------
    cfg : dict
        the brainnets config dictionary
    """
    config.require(cfg,
                   ["all_fnames", "blacklist_fname", "density_range",
                    "node_info_fname", "include_mst", "n_cpus"])
    all_fnames = cfg["all_fnames"]
    arg_list = zip(all_fnames, [cfg] * len(all_fnames))
    ch.run_in_parallel(_link_dist_worker, arg_list, cfg['n_cpus'])


def _link_dist_worker(args):
    fname, cfg = args
    adj_mat, ok_nodes = ch.do_start(fname, cfg['blacklist_fname'])
    all_distances = []
    # load voxeldata
    node_info = dataio.load_mat(cfg['node_info_fname'])['rois'][0]
    node_info = node_info[ok_nodes]
    coords = np.zeros((len(node_info), 3))
    for j in range(len(node_info)):
        coords[j] = np.ravel(node_info['centroidMNI'][j])  # 1 for the coords

    for d in cfg['density_range']:
        net = netgen.make_net_from_unfiltered_data(
            adj_mat, ok_nodes, density=d, include_mst=False, weighted=False)
        es = net.es
        # compute results:
        distances = np.zeros(len(es))
        for j, e in enumerate(es):
            source_coords = coords[e.source]
            target_coords = coords[e.target]
            distances[j] = np.linalg.norm(source_coords - target_coords)
        all_distances.append(distances)
    out_dict = {settings.link_distance_tag: all_distances,
                settings.config_tag: cfg}
    out_fname = fnc.get_ind_fname(fname, cfg, settings.link_distance_tag)
    ch.do_end(fname, out_fname, out_dict)


def comp_paired_common_and_differing_link_distances(cfg):
    """
    For each pair, computes the distances of common and differing
    links.
    """
    config.require(cfg, ["group_1_mat_fnames", "group_2_mat_fnames",
                         "blacklist_fname", "node_info_fname",
                         "density_range", "n_cpus"])
    fnames_group_1 = cfg['group_1_mat_fnames']
    fnames_group_2 = cfg['group_2_mat_fnames']
    assert len(fnames_group_1) == len(fnames_group_2)
    n = len(fnames_group_1)
    arg_list = zip(fnames_group_1, fnames_group_2, [cfg] * n)
    ch.run_in_parallel(_paired_common_and_diff_link_distances_worker,
                       arg_list, cfg['n_cpus'])


def _paired_common_and_diff_link_distances_worker(args):
    fname1, fname2, cfg = args
    adj_mat1, ok_nodes = ch.do_start(fname1, cfg['blacklist_fname'])
    adj_mat2, _ = ch.do_start(fname2, cfg['blacklist_fname'])
    out_dict = {}

    alldistances = []

    for d in cfg['density_range']:
        # print d
        net1 = netgen.make_net_from_unfiltered_data(
            adj_mat1, ok_nodes, density=d, include_mst=True)
        # print "net1"
        net2 = netgen.make_net_from_unfiltered_data(
            adj_mat2, ok_nodes, density=d, include_mst=True)

        net1specific = net1.difference(net2)
        net2specific = net2.difference(net1)
        common2 = net2.difference(net2specific)
        common1 = net1.difference(net1specific)
        for i in range(3):
            assert common1.es[i].source == common2.es[i].source

        distances1 = get_link_distances_for_net(net1specific, cfg)
        distances2 = get_link_distances_for_net(net2specific, cfg)
        distancescommon = get_link_distances_for_net(common1, cfg)
        alldistances.append([distances1, distances2, distancescommon])
    out_dict[settings.link_distance_tag] = alldistances
    out_dict[settings.config_tag] = cfg
    out_fname = fnc.get_paired_fname(fname1, fname2, cfg,
                                     settings.link_distance_common_tag)
    ch.do_end(fname1, out_fname, out_dict)


def comp_consistent_link_distances(cfg):
    """Computes the distances of consistently appearing links
        (amongst the networks defined by the fNames)

        Saves the results to the specified outdata_dir

    Parameters
    ----------
    cfg : dict
        the brainnets config dict
    """
    config.require(cfg, ['all_fnames', 'blacklist_fname', 'node_info_fname',
                         'density_range', 'include_mst'])

    start_mat, ok_nodes = ch.do_start(cfg['all_fnames'][0],
                                      cfg['blacklist_fname'])
    out_dict = {settings.densities_tag: cfg['density_range']}
    link_distances = []
    for d in cfg['density_range']:
        start_net = netgen.make_net_from_unfiltered_data(
            start_mat, ok_nodes, d, include_mst=cfg['include_mst'])
        for i, fName in enumerate(cfg['all_fnames'][1:]):
            mat, ok_nodes = ch.do_start(fName, cfg['blacklist_fname'])
            net = netgen.make_net_from_unfiltered_data(
                mat, ok_nodes, d, include_mst=cfg['include_mst'])
            start_net = start_net.intersection(net)
        link_distances.append(
            get_link_distances_for_net(start_net, cfg)
        )

    out_dict[settings.link_distance_tag] = link_distances

    out_fname = fnc.get_fname(cfg, settings.link_distance_common_tag)
    dataio.save_pickle(out_fname, out_dict)


def get_link_distances_for_net(g, cfg):
    """
    Get link distances for a gwork.

    Parameters
    ----------
    g : igraph.Graph
    cfg : dict
        brainnets config dictionary

    Returns
    -------
    distances : a numpy array
        All distances in the order of the graph's edge sequence.
    """
    node_info = dataio.load_mat(
        cfg['node_info_fname'], squeeze_me=True)["rois"]
    ok_nodes = dataio.get_ok_nodes(cfg['blacklist_fname'])
    coords = node_info['centroidMNI'][ok_nodes]
    distances = np.zeros(len(g.get_edgelist()))
    for j, e in enumerate(g.es):
        source_coords = coords[e.source]
        target_coords = coords[e.target]
        # magic digit two arises from the 2mm NMI brain!
        distances[j] = np.linalg.norm(source_coords - target_coords)
    return distances
