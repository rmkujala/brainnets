from brainnets import compcoms, settings, dataio
import numpy as np
from brainnets import visualizations as viz
from matplotlib import colors
from verkko.plots import colorhelp

import mycfg
cfg = mycfg.get_config()


# preferences to set:
cfg['n_cpus'] = 1
cfg['include_mst'] = False
cfg['density'] = 0.1
# Relatively high density is used as there are not that many nodes (ROI-level networks)
# (in this example, the original network is created on the ROI level to ensure that computations are fast
#  however, we recommend using voxel-level networks in practice.)

# which community detection method to use:
# change e.g. to "infomap" if one wants to use infomap
# other options include methods available in igraph
# if <tag> is used as comdet_tag, igraph should have
# a function called community_<tag>
comdet_tag = settings.louvain_cluster_tag # default
# comdet_tag = "fastgreedy"
# comdet_tag = "label_propagation"
# comdet_tag = "spinglass"
# comdet_tag = "edge_betweenness"
# comdet_tag = "walktrap"
# comdet_tag = "leading_eigenvector" # resulted in some
# Note that you should perhaps also go through the available options
# each of the methods have


n_to_consider = None  # use all partitions by default, other options include e.g. 'best'
if comdet_tag != settings.louvain_cluster_tag:
    assert n_to_consider != "best", "n_to_consider='best' is currently available only for louvain"

# computing prefix for outputting data
ddir = cfg['outdata_dir']
prefix = ddir + "consensus_" + str(cfg['density'])
if cfg['include_mst']:
    prefix = prefix + "_with_mst"
prefix = prefix + "_" + str(comdet_tag)




# compute original communities
print "Computing communities"
# using the Louvain method:

if comdet_tag == settings.louvain_cluster_tag:
    compcoms.comp_louvain_communities(cfg)
else:
    compcoms.comp_communities_igraph(cfg, comdet_tag)


print "Computing consensus communities"
if comdet_tag is None:
    comdet_tag = settings.louvain_cluster_tag
# if n_to_consider is not provided, defaults to using all partitions
consensus_fname0 = (cfg['outdata_dir'] + comdet_tag + "_" + "group1" + "_" +
                    str(cfg['density']) + ".pkl")
consensus_fname1 = (cfg['outdata_dir'] + comdet_tag + "_" + "group2" + "_" +
                    str(cfg['density']) + ".pkl")
compcoms.comp_consensus_partition(cfg, 'group_1_mat_fnames', consensus_fname0,
                                  n_to_consider=n_to_consider, comdet_tag=comdet_tag)
compcoms.comp_consensus_partition(cfg, 'group_2_mat_fnames', consensus_fname1,
                                  n_to_consider=n_to_consider, comdet_tag=comdet_tag)



print "Matching modules algorithmically"
consensus_clu_movie = \
    dataio.load_pickle(consensus_fname0)[comdet_tag]
consensus_clu_rest = \
    dataio.load_pickle(consensus_fname1)[comdet_tag]

consensus_clu_movie, consensus_clu_rest, n_clu_1, n_clu_2 = \
    viz._get_matched_renumbered_communities_for_viz(
        consensus_clu_movie, consensus_clu_rest
    )
uf_consensus_clus = [consensus_clu_movie, consensus_clu_rest]




print "Setting colors for modules"
ok_nodes = dataio.get_ok_nodes(cfg['blacklist_fname'])
filtered_consensus_clus = [
    consensus_clu_movie[ok_nodes], consensus_clu_rest[ok_nodes]]

n_cols = np.maximum(
    np.max(consensus_clu_movie), np.max(consensus_clu_rest)) + 1
# +1 because indexing starts from zero

colors_list = np.array([
    "#e41a1c",  # red
    "#f781bf",  # pink
    "#377eb8",  # blue
    "#95c3dd",  # light blue
    "#1f8f61",  # dark green
    "#a6dc6c",  # light green
    "#984ea3",  # violet
    "#a65628",  # brown
    "#ffff33",  # yellow
    "#ff7f00",  # orange
    "#fbb14e",  # light orange
    "#f88685",  # light red
    "#4daf4a",  # green
])
if n_cols > len(colors_list):
    colors_list = colorhelp.get_arbitrary_n_of_distinct_colors(n_cols)

cc = colors.ColorConverter()
colors_list = cc.to_rgba_array(colors_list)
print colors_list
mod_colors1 = colors_list[np.sort(np.unique(filtered_consensus_clus[0]))]
mod_colors2 = colors_list[np.sort(np.unique(filtered_consensus_clus[1]))]
mod_colors = [mod_colors1, mod_colors2]
both_colors = colors_list


if True:  # plot_alluvial
    print "plotting alluvial"
    fig = viz.plot_alluvial_diagram(
        consensus_clu_rest,  # non-filtered
        consensus_clu_movie,  # non-filtered
        cfg["node_info_fname"],
        cfg['blacklist_fname'],
        ribbon_label_size_lim=15,
        mod_colors1=mod_colors2,
        mod_colors2=mod_colors1,
        mask_large_cluster_w_gray=False
    )
    fig.savefig(prefix + "_alluvial.pdf")

if True:  # plot_on_brain, requires NIFTI toolbox!
    print "plotting consensus modules with slices"
    fig = viz.viz_com_structure_using_slices(
        consensus_clu_movie[ok_nodes], cfg,
        module_colors=mod_colors1
    )
    fig.savefig(prefix + "_movie_on_brain.pdf")

    # filtered coms here:
    fig = viz.viz_com_structure_using_slices(
        consensus_clu_rest[ok_nodes], cfg, module_colors=mod_colors2)
    fig.savefig(prefix + "_rest_on_brain.pdf")


if True:  # plot coarse grained
    from matplotlib import rc
    rc('text', usetex=True)
    # plotting (and computing) coarse-grained networks
    setups = ['movie', 'rest']

    for setup, com, mod_colors_i in zip(setups, uf_consensus_clus, mod_colors):
        fig, cg_mats = viz.comp_and_viz_cluster_diff_matrices(
            com,
            cfg,
            cfg['all_fnames'],
            len(cfg['all_fnames']) / 2,
            vmin=-5,
            vmax=5,
            suptitle="",
            recompute=False,
            mod_colors=mod_colors_i,
            invert_cmap=False,
            return_cg_mats=True,
            titles=setups + ['difference']
        )
        # store results of the coarse-graining:
        out_dict = {"cgmats": cg_mats[::-1],
                    "info"  : "coarse-grained matrices using " \
                              + setup \
                              + " consensus modules, first coarse-grained movie networks," \
                              + " then coarse-grained rest networks"
                    }
        out_fname = cfg["outdata_dir"] + "coarse_grained_matrices_with_" + setup + "_modules"
        dataio.save(out_fname + ".pkl", out_dict)
        dataio.save(out_fname + ".mat", out_dict)
        fig.savefig(prefix + "_" + setup + "_coarse_grained.pdf")