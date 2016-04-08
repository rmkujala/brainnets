import os
import warnings
import numpy as np
from matplotlib.colors import ColorConverter

import settings
CONFIG_KEYS = [
    # file names that are considered when computing networkwise properties
    "all_fnames",
    # name of the blacklist filename
    "blacklist_fname",
    # number of bootstrap samples, required by compstats
    # (and by some plotting routines)
    "bootstrap_samples",
    # coverage of the bootstrap interval, required by compstats
    # (and by some plotting routines)
    "bootstrap_coverage",
    # network density, 0.01 equals to 1\%
    # then the network has density*(n*(n-1))/2 links in it
    "density",
    # range of network densities, if sweeped over a range of network densities
    "density_range",
    # filenames corresponding to the first group,
    # when comparative statistics are made
    "group_1_mat_fnames",
    # filenames corresponding to second group
    # when comparative statistics are made
    "group_2_mat_fnames",
    # group labels, used in plotting
    "group_1_label",
    "group_2_label",
    # group colors, used in plotting
    "group_1_color",
    "group_2_color",
    # include mst, used when thresholding
    "include_mst",
    # which community detection algorithm to use
    # number of iterations used with the Louvain community detection algorithm
    # (Louvain algorithm is stochastic, thus runs can yeild different results)
    "n_it_comdet",
    # number of permutations used in the
    "n_it_permutation",
    # information about the network nodes
    "node_info_fname",
    # where to put all the output data
    "outdata_dir",
    # whether any comparisons are made in a paired fashion or not,
    # used in compstats
    "paired",
    # number of cpus used in the computation
    "n_cpus",
    # node properties
    "node_props",
    # global unweighted properties
    "global_uw_props",
    # global weighted properties
    "global_w_props"
]



def get_default_config():
    """
    Get a default config dictionary, which is filled as much is reasonably
    possible.

    Returns
    -------
    cfg : dict
        default brainnets config dict
    """
    cfg = {}
    return fill_defaults(cfg)


def check_full_config(cfg):
    """
    Checks the validity of all config parameters.
    See :py:func:require.

    Parameters
    ----------
    cfg : dict
        a brainnets config dict to be checked

    Returns
    -------
    valid : bool
        True, if all parameters were ok.

    Raises
    ------
    AssertionError
        If some property is not valid.
    """
    return require(cfg, CONFIG_KEYS)


def fill_defaults(cfg):
    """
    Fills in the default parameters for each key, for which no value is
    currently defined.

    Parameters
    ----------
    cfg : dict
        a brainnets config dict.

    Returns
    -------
    cfg : dict
        the (possibly updated) cfg dict
    """
    for key in CONFIG_KEYS:
        if key not in cfg:
            set_default(cfg, key)
    return cfg


def set_default(cfg, key):
    """
    Fills the default input variables to the extent reasonable
    defaults can be defined.

    Parameters
    ----------

    cfg : dict
        the brainnets config dict
    key : str
        the specific key of the config dictionary to be set as default.

    Returns
    -------
    cfg : dict
        the updated config dict.
    """
    val = ""
    if key is "all_fnames":
        val = None  # no proper default
    elif key is "blacklist_fname":
        val = None  # no proper default
    elif key is "bootstrap_samples":
        val = 10e5  # to be safe
    elif key is "bootstrap_coverage":
        val = 95  # pretty standard
    elif key is "density":
        val = 0.02  # corresponding to 2%, arbitrary
    elif key is "density_range":
        # arbitrary, but reasonable sweep for fine-grained networks:
        val = 2 ** np.linspace(-4, 4, 17) * 0.01
    elif key is "group_1_mat_fnames":
        val = None
    elif key is "group_2_mat_fnames":
        val = None
    elif key is "group_1_label":
        val = "group1"
    elif key is "group_2_label":
        val = "group2"
    elif key is "group_1_color":
        val = "red"
    elif key is "group_2_color":
        val = "blue"
    elif key is "include_mst":
        val = False  # do not include mst by default
    elif key is "n_it_comdet":
        val = 100
    elif key is "paired":
        val = None  # no proper default
    elif key is "outdata_dir":
        try:
            val = os.path.dirname(cfg["all_fnames"][0])
        except:
            try:
                val = os.path.dirname(cfg["group_1_mat_fnames"][0])
            except Exception:
                val = None
    elif key is "n_it_permutation":
        if "paired" in cfg and cfg["paired"] is True:
            val = "all"
        else:
            val = 10e5
    elif key is "node_info_fname":
        val = None  # no proper default
    elif key is "n_cpus":
        val = 1
    elif key is "node_props":
        val = settings.node_props
    elif key is "global_w_props":
        val = settings.global_w_props
    elif key is "global_uw_props":
        val = settings.global_uw_props
    elif key is "com_det_algo":
        val = "louvain"
    cfg[key] = val
    return cfg


def check_property(cfg, key):
    """
    Checks the validity (in a loose sense) of the brainnets
    config dictionary's (``cfg``) value corresponding to ``key``.key

    Parameters
    ----------
    cfg : dict
        brainnets config dict
    key : str
        the key of the dictionary, which should be checked

    Returns
    -------
    valid : bool
        True, if all parameters were ok.

    Raises
    ------
    AssertionError
        If the value of ``cfg[key]`` was not valid.
    """
    key_in_config = (key in cfg) and cfg[key] is not None
    assert key_in_config, "config key \"" + str(key) + "\" is not specified"

    props_tags = ["node_props", "global_uw_props", "global_w_props"]
    settings_props_list = [settings.node_props,
                           settings.global_uw_props,
                           settings.global_w_props]

    def _asm():
        return "config key " + str(key) + " is not properly set," + \
            "\n value: " + str(val)

    val = cfg[key]
    # check the actual assertions
    if key is "density":
        assert (0 <= val and val <= 1), _asm()

    elif key in ["group_1_color", "group_2_color"]:
        try:
            cc = ColorConverter()
            cc.to_rgba(val)
        except ValueError:
            raise AssertionError(_asm())

    elif key in ["all_fnames", "group_1_mat_fnames", "group_2_mat_fnames",
                 "node_info_fname", "blacklist_fname", "outdata_dir"]:
        if type(val) is not list:
            val = [val]
        for fname in val:
            assert os.path.exists(fname), "check key: " + key
    elif key == "bootstrap_samples":
        try:
            float(val)
            assert val > 10, _asm() + " (may be too low)"
        except:
            assert val == "all", _asm()
    elif key is "bootstrap_coverage":
        assert val > 0 and val < 101, _asm()
    elif key is "density_range":
        val = np.array(val)
        assert ((val <= 1) * (val >= 0)).all(), _asm()
    elif key in ["group_1_label", "group_2_label"]:
        assert len(val) > 0 and type(val) is str, _asm()
    elif key in ["paired", "include_mst"]:
        assert type(val) is bool, _asm()
    elif key is "n_it_comdet":
        assert val >= 1, _asm()
    elif key is "n_it_permutation":
        if type(val) is str:
            assert val == "all", _asm()
        else:
            assert val >= 4
            if key is "n_it_permutation" and val < 100:
                warnings.warn("Low number of permutations for key " + str(key)
                              + ": " + str(val))
    elif key is "n_cpus":
        assert (type(val) is int) and val >= 1
    elif key in props_tags:
        i = props_tags.index(key)
        props_should_be_in = settings_props_list[i]
        for v in val:
            assert v in props_should_be_in, _asm()
    elif key not in CONFIG_KEYS:
        warnings.warn('Omitting the check of custom config param ' + str(key))
    else:
        warnings.warn("The validity of " + str(key) +
                      " was not properly checked ")
    return True


def require(cfg, keys):
    """
    Checks the validity of config dictionary's values corresponding to the
    ``keys``.

    Parameters
    ----------
    cfg : dict
        the brainnets config dictionary
    keys : list of strings
        the validity of the keys to be checked


    Returns
    -------
    valid : bool
        True, if all parameters were ok.

    Raises
    ------
    AssertionError
        If at least one of the keys is not properly defined.
    """
    for key in keys:
        check_property(cfg, key)
    return True
