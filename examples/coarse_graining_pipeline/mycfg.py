from brainnets import config, settings
import os

"""
See config (and the compprops etc. functions) for further explanations of the
the meanings config parameters
"""

# change this to the correct path
# settings.path_to_NIFTI_toolbox = "/path/to/NIFTI"


def get_config():
    path_to_this_dir = os.path.dirname(os.path.realpath(__file__))
    input_dir = path_to_this_dir + "/input_data/"

    movie_fnames = [input_dir + "movie_HO_" + str(i) + ".mat" for i in range(1, 14)]
    rest_fnames = [input_dir + "rest_HO_" + str(i) + ".mat" for i in range(1, 14)]

    # the values here should not be considered as the default ones
    # (defaults can be obtained from config.py:get_default_config()
    cfg = {
        "all_fnames": movie_fnames + rest_fnames,
        "blacklist_fname": input_dir + "blacklist_bool.mat",
        "bootstrap_samples": 1000,
        "bootstrap_coverage": 95,
        "density": 0.02,
        "density_range": [0.03, 0.1, 0.5],
        "group_1_mat_fnames": movie_fnames,
        "group_2_mat_fnames": rest_fnames,
        "group_1_label": "movie",
        "group_2_label": "rest",
        "group_1_color": "red",
        "group_2_color": "blue",
        "include_mst": True,  # used in thresholding
        "n_it_comdet": 2,
        "n_it_permutation": "all",
        "paired": True,
        "n_cpus": 1,
        "node_info_fname": input_dir + "HO_2mm_rois.mat",
        "outdata_dir": path_to_this_dir + "/output/",
        "node_props": settings.node_props,
        "global_uw_props": settings.global_uw_props,
        "global_w_props": settings.global_w_props,
        "n_nodes": 30
    }
    config.check_full_config(cfg)
    return cfg
