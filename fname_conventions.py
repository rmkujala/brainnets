# from brainnets import settings
import settings


def extract_basename(fname):
    """
    Extracts the basename of a filename, which can be used as a prefix for
    resulting data files. Allows also for "." in the basename.

    Parameters
    ----------
    fname : str

    Returns
    -------
    basename : str


    Examples
    --------

    >>> fname = "foo/bar/foo.bar.tadaa"
    >>> extract_basename(fname)
    'foo.bar'

    """
    splitted = fname.split("/")[-1].split(".")[:-1]
    basename = ""
    for i, s in enumerate(splitted):
        if i >= 1:
            basename += "."
        basename += s
    return basename


def get_ind_fname(fname, cfg, prop_prefix):
    postfix = ""
    if prop_prefix in ["node_props", 'louvain_clusters']:
        postfix += "_density_" + str(cfg['density'])
        if cfg['include_mst']:
            postfix = postfix + "_with_mst"
    return cfg['outdata_dir'] + prop_prefix + "_" + extract_basename(fname) \
        + postfix + ".pkl"


def get_group_fname(cfg, prop_prefix, group_number):
    """
    Get a filename for results conserning one group.

    Parameters
    ----------
    cfg : dict
        the brainnets config dictionary
    prop_prefix : str
        one of the properties listed in brainnets.settings
    group_number : int
        the number of the group (starting from zero)

    Returns
    -------
    fname : str
        name of the file where to write the results
    """
    if prop_prefix in [settings.louvain_consensus_tag,
                       settings.scaled_inclusivity_tag,
                       settings.louvain_consensus_si_tag]:
        postfix = "group_" + \
            str(int(group_number)) + "_density_" + str(cfg['density'])
    else:
        raise Exception("group fname undefined for " + prop_prefix)
    if cfg['include_mst']:
        postfix = postfix + "_with_mst"
    return cfg['outdata_dir'] + prop_prefix + "_" + postfix + ".pkl"


def get_fname(cfg, prop_prefix):
    """
    Get a filename for results considering multiple networks.
    """
    # file names for multiple networks
    postfix = prop_prefix
    if prop_prefix is settings.louvain_cluster_tag:
        postfix = postfix + "_" + str(cfg['density'])
    elif prop_prefix is settings.louvain_consensus_tag:
        postfix = postfix + "_" + str(cfg['density'])
    if cfg['include_mst']:
        postfix = postfix + "_with_mst"
    return cfg['outdata_dir'] + postfix + ".pkl"


def get_paired_fname(fname1, fname2, cfg, prop_prefix):
    """
    Get filename for paired comparisons.
    """
    postfix = ""
    if cfg['include_mst']:
        postfix = "_with_mst"
    return cfg['outdata_dir'] + prop_prefix + extract_basename(fname1) + \
        "_vs_" + extract_basename(fname2) + postfix + ".pkl"


def get_stats_fname(cfg, prop_prefix):
    """
    Get stats fname for some property

    Parameters
    ----------
    cfg : dict
        the brainnets config dictionary
    prop_prefix : str
        one of the properties listed in brainnets.settings

    Returns
    -------
    fname : str
        the stats filename
    """
    return cfg['outdata_dir'] + prop_prefix + "_stats.pkl"
