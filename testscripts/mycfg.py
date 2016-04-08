from brainnets import config, settings
import os

"""
See config (and the compprops etc. functions) for further explanations of the
the meanings config parameters
"""


# def get_old_config():

#     brainnets_dir = os.path.dirname(config.__file__)
#     ddir = brainnets_dir + "/testdata/"

#     movie_fnames = [
#         ddir + "movie_small_" + str(i) + ".mat" for i in range(1, 14)]
#     rest_fnames = [ddir + "rest_small_" +
#                    str(i) + ".mat" for i in range(1, 14)]

# the values here should not be considered as the default ones
# (defaults can be obtained from config.py:get_default_config()
#     cfg = {
#         "all_fnames": movie_fnames + rest_fnames,
#         "blacklist_fname": ddir + "jointblacklist_small.mat",
#         "bootstrap_samples": 100,
#         "bootstrap_coverage": 50,
#         "density": 0.05,
#         "density_range": [0.1, 0.03, 0.5],
#         "group_1_mat_fnames": movie_fnames,
#         "group_2_mat_fnames": rest_fnames,
#         "group_1_label": "movie",
#         "group_2_label": "rest",
#         "group_1_color": "blue",
#         "group_2_color": "red",
# "include_mst": False,  # False,  # used in thresholding
#         "n_it_comdet": 10,
#         "n_it_permutation": "all",
#         "paired": True,
#         "n_cpus": 1,
#         "node_info_fname": ddir + "rois6mm_final_plus_IF2.mat",
#         "outdata_dir": ddir + "output/",
#         "node_props": settings.node_props,
#         "global_uw_props": settings.global_uw_props,
#         "global_w_props": settings.global_w_props,
#         "n_nodes": 30
#     }
#     return cfg


def get_config():

    brainnets_dir = os.path.dirname(config.__file__)
    ddir = brainnets_dir + "/testdata/"

    movie_fnames = [ddir + "movie_HO_" + str(i) + ".mat" for i in range(1, 14)]
    rest_fnames = [ddir + "rest_HO_" + str(i) + ".mat" for i in range(1, 14)]

    # the values here should not be considered as the default ones
    # (defaults can be obtained from config.py:get_default_config()
    cfg = {
        "all_fnames": movie_fnames + rest_fnames,
        "blacklist_fname": ddir + "blacklist_bool.mat",
        "bootstrap_samples": 100,
        "bootstrap_coverage": 50,
        "density": 0.5,
        "density_range": [0.03, 0.1, 0.5],
        "group_1_mat_fnames": movie_fnames,
        "group_2_mat_fnames": rest_fnames,
        "group_1_label": "movie",
        "group_2_label": "rest",
        "group_1_color": "blue",
        "group_2_color": "red",
        "include_mst": False,  # False,  # used in thresholding
        "n_it_comdet": 10,
        "n_it_permutation": "all",
        "paired": True,
        "n_cpus": 1,
        "node_info_fname": ddir + "HO_2mm_rois.mat",
        "outdata_dir": ddir + "output/",
        "node_props": settings.node_props,
        "global_uw_props": settings.global_uw_props,
        "global_w_props": settings.global_w_props,
        "n_nodes": 30
    }

    return cfg
