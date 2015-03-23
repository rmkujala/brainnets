# brainnets imports
import dataio
import settings
import gencomps
import netgen
import config
import fname_conventions as fnc
import comp_helpers as ch


"""
| The structure of the module could be improved, although it works.
| Everything is made so that the parallellism part
| was needed to code only once.
| Functions starting with an underscore are (at least originally) meant to be
| used only internally.
"""

# Node properties computations:


def comp_node_props(cfg):
    """
    Computes node properties for ``cfg["all_fnames"]`` for a given network
    density

    Parameters
    ----------
    cfg : dict
        brainnets config dictionary
    """
    config.require(
        cfg,
        ["all_fnames", "blacklist_fname", "density", "include_mst", "n_cpus",
         "node_props"]
    )
    all_fnames = cfg["all_fnames"]
    n = len(all_fnames)
    arg_list = zip(all_fnames, [cfg] * n)
    ch.run_in_parallel(_node_prop_worker, arg_list, cfg["n_cpus"])


def _node_prop_worker(args):
    fname, cfg = args
    adj_mat, ok_nodes = ch.do_start(fname, cfg["blacklist_fname"])
    results = gencomps.get_node_props_from_mat(
        adj_mat, ok_nodes, cfg["density"],
        props=cfg["node_props"],
        include_mst=cfg["include_mst"])
    results[settings.config_tag] = cfg
    out_fname = fnc.get_ind_fname(fname, cfg, "node_props")
    ch.do_end(fname, out_fname, results)


def compute_global_properties(cfg, weighted=False):
    """
    Computes global network properties for ``cfg["all_fnames"]`` for a given
    network density

    Parameters
    ----------
    cfg : dict
        brainnets config dictionary
    weighted : bool
        if True, the network is considered as weighted.
        if False, network is unweighted
    """
    config.require(
        cfg,
        ["all_fnames", "blacklist_fname", "density_range", "include_mst",
         "n_cpus",
         "global_w_props"]
    )
    n = len(cfg['all_fnames'])
    arg_list = zip(cfg['all_fnames'], [cfg] * n)
    if weighted:
        worker = _global_w_prop_worker
    else:
        worker = _global_uw_prop_worker
    ch.run_in_parallel(worker, arg_list, cfg['n_cpus'])


def _global_w_prop_worker(args):
    fname, cfg = args
    adj_mat, ok_nodes = ch.do_start(fname, cfg['blacklist_fname'])
    results = gencomps.get_global_props_for_density_range(
        adj_mat,
        ok_nodes,
        cfg['density_range'],
        cfg['global_w_props'],
        True,  # weighted
        cfg['include_mst'],
    )
    results[settings.config_tag] = cfg
    out_fname = fnc.get_ind_fname(fname, cfg, "global_w_props")
    ch.do_end(fname, out_fname, results)


def _global_uw_prop_worker(args):
    fname, cfg = args
    adj_mat, ok_nodes = ch.do_start(fname, cfg['blacklist_fname'])
    results = gencomps.get_global_props_for_density_range(
        adj_mat,
        ok_nodes,
        cfg['density_range'],
        cfg['global_uw_props'],
        False,  # weighted
        cfg['include_mst'],
    )
    results[settings.config_tag] = cfg
    out_fname = fnc.get_ind_fname(fname, cfg, "global_uw_props")
    ch.do_end(fname, out_fname, results)


#
# LINK SIMILARITIES
#


def _comp_link_sim_mat_worker(args):
    density, cfg = args
    nets = []
    for fname in cfg['all_fnames']:
        adj_mat, ok_nodes = ch.do_start(fname, cfg['blacklist_fname'])
        nets.append(netgen.make_net_from_unfiltered_data(
            adj_mat, ok_nodes, density, weighted=False,
            include_mst=cfg['include_mst']))
    return gencomps.comp_link_sim_mat(nets)


def comp_link_sim_matrices(cfg):
    """
    Compute link similarity matrices for given network densities.
    specified in ``cfg['density_range']`` across
    ``cfg['all_fnames']``
    Results are saved in ``cfg['outdata_dir'])``

    Parameters
    ----------
    cfg : dict
        brainnets config dict
    """
    config.require(cfg,
                   ["all_fnames", "blacklist_fname", "density_range",
                    "include_mst", "n_cpus"])
    densities = cfg['density_range']
    arg_list = zip(densities, [cfg] * len(densities))
    link_sim_mat_list = ch.run_in_parallel(
        _comp_link_sim_mat_worker, arg_list,
        cfg["n_cpus"], chunksize=1)
    out_dict = {settings.densities_tag: densities,
                settings.common_links_tag: link_sim_mat_list,
                settings.config_tag: cfg}
    dataio.save_pickle(fnc.get_fname(cfg, settings.common_links_tag), out_dict)
