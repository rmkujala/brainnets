import os
import hashlib
# third party, numpy and matplotlib
import numpy as np
import nibabel as nib
from matplotlib import cm, colors
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
# brainnets
from brainnets import settings, aux, dataio
from brainnets import netgen, communities
from brainnets import fname_conventions as fnc
# verkko
from verkko.plots import alluvial
from verkko.plots import colorhelp
from verkko.permtests import measures
from verkko.permtests import ptests
from verkko.permtests import corrections


def viz_data_on_brain_using_slices(
        in_data_nii_fname="",
        T=1E-300,
        slices="xyz",
        bg_nii_fname="/usr/share/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz",
        colormap=None,
        slice_step=1,
        cmap_lims=None):
    """
    Plotting the indata on top of the MNI template brain.

    Parameters
    ----------
    in_data_nii_fname : str
        path to the .nii image which will be visualized
    T : float
        All values with abs(val) < T are not shown.
    slices : str
        string of all of the slice directions, e.g. ("yz") for y and z slices
    bg_nii_fname : str
        which background image to use
    colormap: matplotlib colormap
        Colormap to use, if None, a suitable colormap for
        visualizing positiv and negative t-values is used.
    cmap_lims : tuple (or list)
        cmap_lims[0] and cmap_lims[1] correspond to the the minimum
        and maximum values of the colorbar
    """
    in_data = nib.load(in_data_nii_fname)
    in_vol = in_data.get_data()
    bg_data = nib.load(bg_nii_fname)
    bg_vol = bg_data.get_data()
    in_data_min = np.min(in_vol)
    in_data_max = np.max(in_vol)
    in_data_abs_max = np.maximum(np.abs(in_data_min), in_data_max)

    using_default_colormap = False
    if colormap is None:
        using_default_colormap = True
        lowcut = (0.5 - T / in_data_abs_max * 0.5)
        highcut = (0.5 + T / in_data_abs_max * 0.5)

        # colormap combination of builtin winter  & autumn colormaps
        cdict = {
            'blue': ((0.0, 1.0, 1.0),
                     (lowcut, 0.5, 0.7),
                     (highcut, 0.7, 0.0),
                     (1.0, 0.0, 0.0)),
            'green': ((0.0, 0.0, 0.0),
                      (lowcut, 1.0, 0.7),
                      (highcut, 0.7, 0.0),
                      (1.0, 1.0, 1.0)),
            'red': ((0.0, 0.0, 0.0),
                    (lowcut, 0.0, 0.7),
                    (highcut, 0.7, 1.0),
                    (1.0, 1.0, 1.0)),
            'alpha': ((0.0, 0.0, 0.0),
                      (lowcut, 0.0, 0.0),
                      (highcut, 0.0, 1.0),
                      (1.0, 1.0, 0.0))
        }
        colormap = LinearSegmentedColormap('custom', cdict, N=1024)

    step_size = slice_step * 3
    slice_ranges = {'x': range(16, 78, step_size),  # 84,step_size),
                    'y': range(12, 95, step_size),
                    'z': range(8, 76, step_size)
                    }
    slicedir2dim = {'x': 0, 'y': 1, 'z': 2}

    in_vol = np.array(in_vol)

    if cmap_lims is not None:
        vmin = cmap_lims[0]
        vmax = cmap_lims[1]
    else:
        vmax = in_data_abs_max
        vmin = -in_data_abs_max

    for slice_dir in slices:
        slice_range = slice_ranges[slice_dir]
        n_slices = len(slice_range)
        in_vol_view = in_vol.swapaxes(0, slicedir2dim[slice_dir])
        bg_vol_view = bg_vol.swapaxes(0, slicedir2dim[slice_dir])
        nrows = n_slices / 7 + bool(n_slices % 7)
        ncols = n_slices / nrows + bool(n_slices % nrows) + 1
        gs = gridspec.GridSpec(nrows, ncols)
        fig = plt.figure(figsize=(15, 6))
        #    12, step=3
        for j, i in enumerate(slice_range):  # 0,len(in_vol[0])):#,10):
            # print j%(ncols-1), j/(ncols-1)
            ax = fig.add_subplot(gs[j / (ncols - 1), j % (ncols - 1)])

            # axes have been swapped accordingly before
            Z = in_vol_view[i, :, :]

            if slice_dir in "yz":
                bg = bg_vol_view[i, :, :]
            else:
                bg = bg_vol_view[-i, :, :]
            cond = (-T < Z) & (Z < T)
            Z_mask = np.ma.masked_where(cond, Z)
#        else:
#                cond = (Z < 0)
#                Z_mask= np.ma.masked_where(cond, Z)

            if slice_dir in "x":
                ax.imshow(
                    bg.T, cmap="gray", origin="lower", interpolation="nearest")
                im = ax.imshow(Z_mask.T, cmap=colormap, origin="lower",
                               interpolation='nearest', alpha=1.0,
                               vmax=vmax, vmin=vmin)
            if slice_dir in "y":
                ax.imshow(
                    bg.T[:, ::-1], cmap="gray", origin="lower",
                    interpolation="nearest")
                im = ax.imshow(Z_mask.T, cmap=colormap, origin="lower",
                               interpolation='nearest', alpha=1.0,
                               vmax=vmax, vmin=vmin)
            if slice_dir == "z":
                ax.imshow(
                    bg[:, ::-1], cmap="gray", origin="lower",
                    interpolation="nearest")
                im = ax.imshow(Z_mask, cmap=colormap, origin="lower",
                               interpolation="nearest", vmax=vmax, vmin=vmin)

            # add left/right texts
            if slice_dir == 'x':
                if j == 0:
                    ax.text(0.1, 0.5, 'L', horizontalalignment='center',
                            verticalalignment='center', transform=ax.transAxes,
                            color="w", weight="bold")
                if j == len(slice_range) - 1:
                    ax.text(0.9, 0.5, 'R', horizontalalignment='center',
                            verticalalignment='center', transform=ax.transAxes,
                            color="w", weight="bold")
            if slice_dir in 'yz':
                ax.text(0.1, 0.9, 'L', horizontalalignment='center',
                        verticalalignment='center', transform=ax.transAxes,
                        color="0.5", weight="bold")
                ax.text(0.9, 0.9, 'R', horizontalalignment='center',
                        verticalalignment='center', transform=ax.transAxes,
                        color="0.5", weight="bold")

            dim = slicedir2dim[slice_dir]
            slice_coords = np.zeros(3, dtype=np.int16)
            slice_coords[dim] = i
            mni_coords = aux.space2MNI(slice_coords)[dim]
            ax.text(0.97, 0.1, str(mni_coords), horizontalalignment='right',
                    verticalalignment='center', transform=ax.transAxes,
                    color="0.5", weight="bold")

            ax.set_xticks([])
            ax.set_yticks([])

        # colorbar axes
        cax = fig.add_axes([0.9, 0.1125, 0.02, 0.775])
        cbar = plt.colorbar(im, cax=cax)
        if not using_default_colormap:
            maxZ = np.nanmax(in_vol_view)
            print maxZ
            if maxZ - int(maxZ) < 1e-6:
                print maxZ
                if maxZ < 16:
                    ticks = np.arange(1, maxZ + 2)
                    cbar.set_ticks(ticks)
        plt.tight_layout(w_pad=0.2)
    return fig


def viz_nodewise_vals_using_slices(
        node_val_array, node_info_fname, blacklist_fname,
        vals_blacklist_filtered=False,
        default_value=0, slicedirs="xyz", cmap=None,
        cmap_lims=None, T=1e-300, slice_step=1):
    """
    Plots the node values (node_val_array) using different brain slices.

    Parameters
    ----------
    node_val_array : numpy array
        1D numpy array containing the values for each node,
        by default the values should not be blacklist filtered
    vals_blacklist_filtered : bool
        True if node_val_array is blacklist filtered
    blacklist_fname : str
        path to the ".mat" file containing the node blacklist
    node_info_fname : str
        name of the ".mat" file containing the node information
    default_value : node_val_array.dtype
        default value for non-blacklisted nodes
    cmap : a matplotlib colormap
        colormap to be used in the visualization
    cmap_lims : tuple of floats
        the min and max of the colormap range
    T : float
        threshold value used in the visualization
        (to omit exactly zero values, the code needs improvement)
    slice_step : int
        how many 2mm slices are between images
    """
    if vals_blacklist_filtered:
        ok_nodes = dataio.get_ok_nodes(blacklist_fname)
        node_val_array = dataio.expand_1D_node_vals_to_non_blacklisted_array(
            node_val_array, ok_nodes, default_value=default_value)
    pid_str = str(os.getpid())
    nii_fname = "/tmp/" + pid_str + ".nii"

    dataio.out_put_node_prop_to_nii(
        node_val_array,
        "/tmp/" + pid_str + ".nii",
        node_info_fname,
        blacklist_fname)
    fig = viz_data_on_brain_using_slices(
        nii_fname, slices=slicedirs,
        colormap=cmap, cmap_lims=cmap_lims, T=T, slice_step=slice_step)
    os.remove(nii_fname)
    return fig

# works


def viz_node_props_for_ind(fname, prop_tag, cfg, savefig=True):
    """
    Visualizes node properties for one subject.

    Parameters
    ----------
    fname : str
        the original mat file containing the correlation matrix
    prop_tag : str
        the tag of the property
    cfg : dict
        brainnets config dictionary

    Returns
    -------
    fig : a matplotlib.Figure
        the figure with the brain properties plotted on the brain
    """
    data = dataio.load_pickle(fnc.get_ind_fname(fname, cfg, "node_props"))
    node_val_array = data[prop_tag]
    fig = viz_nodewise_vals_using_slices(
        node_val_array, cfg['node_info_fname'],
        cfg['blacklist_fname'], default_value=0
    )

    if savefig:
        basename = (fnc.extract_basename(fname) + "_" + str(prop_tag) + "_"
                    + str(cfg["density"]))
        fig.savefig(cfg['outdata_dir'] + basename + ".pdf", format="pdf")
    return fig


def viz_com_stru_for_ind(fname, cfg, i=0, savefig=True):
    """
    Visualizes the community structure for one subject.

    Parameters
    ----------
    fname : str
        the original mat file containing the correlation matrix
    cfg : dict
        brainnets config dictionary
    i : int
        which com structure to plot if multiple communities have been computed

    Returns
    -------
    fig : matplotlib.Figure
        the figure with the brain properties plotted on the brain
    """
    data = dataio.load_pickle(
        fnc.get_ind_fname(fname, cfg, settings.louvain_cluster_tag)
    )
    com = data[settings.louvain_cluster_tag][i]
    ok_nodes = dataio.get_ok_nodes(cfg['blacklist_fname'])
    com = com[ok_nodes]
    fig = viz_com_structure_using_slices(
        com, cfg, slicedirs="x",
        mask_large_cluster_w_gray=True
    )
    if savefig:
        basename = (fnc.extract_basename(fname) + "_" +
                    str(settings.louvain_cluster_tag) + "_" +
                    str(cfg["density"]) + "_" + str(i))
        fig.savefig(cfg['outdata_dir'] + basename + ".pdf", format="pdf")
    return fig


def viz_com_structure_from_file(com_fname, cfg, module_colors=None, **kwargs):
    """
    A convenience function for

    Parameters
    ----------
    com_fname : str
        path to the filename which is going to be visualized
    cfg : dict
        brainnets config dict
    kwargs : optional, see :py:func:`viz_com_structure_using_slices`
    """
    com = dataio.load(com_fname)[settings.louvain_cluster_tag]
    ok_nodes = dataio.get_ok_nodes(cfg['blacklist_fname'])
    com = com[ok_nodes]
    return viz_com_structure_using_slices(com, cfg, module_colors=module_colors, **kwargs)


def viz_com_structure_using_slices(filtered_comstructure, cfg, module_colors=None, slicedirs="x",
                                   mask_large_cluster_w_gray=False):
    """
    Vizualizes a given unfiltered community structure in brain slices.

    Parameters
    ----------
    filtered_comstructure : list / numpy array
        membership list, contains values only for the ok_nodes
    cfg : dict
        brainnets confic dict
    slicedirs : str
        e.g. "yz" means that slices whose tangent is aligned with "y" and
        "z" axes are shown
    mask_large_cluster_with_gray : bool
        whether to mask one big cluster (>=0.5*N) with gray

    Returns
    -------
    fig : matplotlib figure

    See also
    --------
    viz_com_structure_from_file : give a filename to a com structure
        instead of a filtered_comstructure
    """
    filtered_comstructure = np.array(filtered_comstructure, dtype=np.int32)
    sorted_clu_labels = np.setdiff1d(
        filtered_comstructure, [settings.undef_clu_label])

    newComStructure = np.zeros(filtered_comstructure.shape)
    for i, label in enumerate(sorted_clu_labels):
        indices = (filtered_comstructure == label)
        newComStructure[indices] = i

    n_colors = len(sorted_clu_labels)
    if module_colors is None:
        cols = get_module_colors(np.sort(np.unique(sorted_clu_labels)))
    else:
        cols = module_colors
    # newComStructure = communities.make_zero_indexed_clustering(
    #     filtered_comstructure)

    # mask some cluster with gray?:
    if mask_large_cluster_w_gray:
        tmpComStructure = newComStructure[
            newComStructure != settings.undef_clu_label]
        clusizes = communities.get_module_sizes(np.array(tmpComStructure))
        if mask_large_cluster_w_gray:
            mask_large_cluster_with_gray(cols, clusizes)

    listed_color_map = colors.ListedColormap(cols, N=len(cols))
    listed_color_map.set_bad(color="k", alpha=0.0)

    fig = viz_nodewise_vals_using_slices(
        newComStructure + 1, cfg['node_info_fname'], cfg['blacklist_fname'],
        slicedirs=slicedirs, vals_blacklist_filtered=True, default_value=0,
        cmap=listed_color_map,
        cmap_lims=[0.5, n_colors + 0.5])
    return fig

# PLOTTING THE BRAIN NODE PROPERTIES USING BRAIN SLICES


def _get_matched_renumbered_communities_for_viz(
        clus_1,
        clus_2,
        undef_clu_label=settings.undef_clu_label):
    """
    Computes new matched and renumbered (starging from zero)
    communities.

    Parameters
    ----------
    clus_1 : np.array
        membership list, with non-ok nodes marked with
        :py:attr:`settings.undef_clu_label
    clus_2 : np.array
        membership list, with non-ok nodes marked with
        :py:attr:`settings.undef_clu_label

    Returns
    -------
    filtered_clus_1 : np.array
    filtered_clus_2 : np.array
    n_clu_1 : int
        number of clusters in the first partition
    n_clu_2 : int
        number of clusters in the second partition
    """
    clus_1 = np.array(clus_1, dtype=np.int32)
    clus_2 = np.array(clus_2, dtype=np.int32)
    ok_nodesMovie = (clus_1 != undef_clu_label)
    ok_nodesRest = (clus_1 != undef_clu_label)
    assert (ok_nodesMovie == ok_nodesRest).all()
    ok_nodes = ok_nodesMovie

    filtered_clus_1 = clus_1[clus_1 != undef_clu_label] + 1
    filtered_clus_2 = clus_2[clus_2 != undef_clu_label] + 1

    n_clu_1 = len(np.unique(filtered_clus_1))
    n_clu_2 = len(np.unique(filtered_clus_2))

#    if n_clu_1 > n_clu_2:
#        filtered_clus_2 = \
#            gencomps.matchClustersHungarianAlgo(filtered_clus_1,
#                                                filtered_clus_2)
#    else:
#        filtered_clus_1 = \
#             gencomps.matchClustersHungarianAlgo(filtered_clus_2,
#                                                 filtered_clus_1)
# in place of the above(?)

    filtered_clus_1, filtered_clus_2 = \
        communities.match_clusters_greedy(filtered_clus_1,
                                          filtered_clus_2
                                          )

    all_clu_nums = np.union1d(filtered_clus_1, filtered_clus_2)

    # this could maybe be simplified
    filtered_clus = [filtered_clus_1, filtered_clus_2]
    for i in range(len(filtered_clus)):
        clustering = filtered_clus[i]
        new_clu = np.ones(len(clustering), dtype=np.int64) * -1000
        for j, clu_num in enumerate(all_clu_nums):
            clunum_indices = (clustering == clu_num)
            new_clu[clunum_indices] = j
        # expands only if ok_nodes is not all
        filtered_clus[i] = dataio.expand_1D_node_vals_to_non_blacklisted_array(
            new_clu, ok_nodes, default_value=undef_clu_label)

    filtered_clus_1 = filtered_clus[0]
    filtered_clus_2 = filtered_clus[1]
    return filtered_clus_1, filtered_clus_2, n_clu_1, n_clu_2


def plot_alluvial_diagram(consensus_clu_movie,
                          consensus_clu_rest,
                          node_info_fname,
                          blacklist_fname,
                          ribbon_label_size_lim=15,
                          stab_mask_1=None,
                          stab_mask_2=None,
                          mask_large_cluster_w_gray=True,
                          mod_colors1=None,
                          mod_colors2=None):
    """
    Plots the alluvial diagram between two community structures
    (com_struct_fname_1 & com_struct_fname_2), see what is expected from
    the files in the code below.

    Parameters
    ----------
    com_struct_fname_1 : str
        a .pkl file containing the first non-blacklisted community structure
    com_struct_fname_2 : str
        a .pkl file containing the second non-blacklisted community structure
    node_info_fname : str
        path to the node info fname (.mat)
    blacklist_fname :
        path to the blacklist
    ribbon_label_size_lim : int
        how many labels of some type ("label") is required for adding a label
    stab_mask_1 : array of bools
        a (non-blacklisted) cluster mask for the stable nodes in
        com_struct_fname_1. If i is stable, then stab_mask_1[i] = True (?)
    stab_mask_2 : array of bools
        a (non-blacklisted) cluster mask for the stable nodes in
        com_struct_fname_1. If i is stable, then stab_mask_1[i] = True (?)
    mask_large_cluster_with_gray : bool
        if cluster's size is more than half of the total number of nodes
        then it's color is forced to be grey

    Returns
    -------
    fig : the matplotlib figure
    """
    if stab_mask_1 is not None and stab_mask_2 is not None:
        plot_stabilities = True
    else:
        plot_stabilities = False

    rois = dataio.load_mat(node_info_fname, squeeze_me=True)["rois"]
    ok_nodes = dataio.get_ok_nodes(blacklist_fname)
    okRois = rois[ok_nodes]
    # okRoiAalLabels = okRois['better_label']
    okRoiAalLabels = okRois['abbr']

    if mask_large_cluster_with_gray:
        max_val = np.max([consensus_clu_movie, consensus_clu_rest])
        for conclu in [consensus_clu_movie, consensus_clu_rest]:
            clu_labels = np.setdiff1d(conclu, [settings.undef_clu_label])
            for i, item in enumerate(clu_labels):
                mod_nodes = (item == conclu)
                if np.sum(mod_nodes) > 0.5 * len(mod_nodes) \
                        and item != max_val:
                    conclu[mod_nodes] = max_val + 1

    consensus_clu_movie = consensus_clu_movie[ok_nodes]
    consensus_clu_rest = consensus_clu_rest[ok_nodes]

    n_clu_1 = len(np.unique(consensus_clu_movie))
    n_clu_2 = len(np.unique(consensus_clu_rest))

    if plot_stabilities:
        stab_mask_1 = stab_mask_1[ok_nodes]
        stab_mask_2 = stab_mask_2[ok_nodes]
        stableRibbonSizes1 = np.zeros((n_clu_1, n_clu_2))
        stableRibbonSizes2 = np.zeros((n_clu_1, n_clu_2))
    else:
        stab_mask_1 = None
        stab_mask_2 = None
        stableRibbonSizes1 = None
        stableRibbonSizes2 = None

    if mod_colors1 is None or mod_colors2 is None:
        mod_colors1 = get_module_colors(
            np.sort(np.unique(consensus_clu_movie)))
        mod_colors2 = get_module_colors(np.sort(np.unique(consensus_clu_rest)))

    # compute modulesizes & ribbon_size_mat
    moduleSizes1 = np.zeros(n_clu_1)
    moduleSizes2 = np.zeros(n_clu_2)

    np.zeros((n_clu_1, n_clu_2))
    ribbon_label_mat = np.zeros((n_clu_1, n_clu_2), dtype=object)
    ribbon_size_mat = np.zeros((n_clu_1, n_clu_2), dtype=float)

    clu_labels = np.setdiff1d(consensus_clu_movie, [settings.undef_clu_label])
    for i, item in enumerate(clu_labels):
        moduleNodes1 = (item == consensus_clu_movie)
        moduleSizes1[i] = np.sum(item == consensus_clu_movie)

        for j, jtem in enumerate(np.sort(np.unique(consensus_clu_rest))):
            moduleNodes2 = (jtem == consensus_clu_rest)
            moduleSizes2[j] = np.sum(moduleNodes2)
            commonNodes = moduleNodes2 * moduleNodes1
            nCommonNodes = np.sum(moduleNodes2 * moduleNodes1)
            ribbon_size_mat[i, j] = nCommonNodes

            if plot_stabilities:
                stableRibbonSizes1[i, j] = np.sum(commonNodes * stab_mask_1)
                stableRibbonSizes2[i, j] = np.sum(commonNodes * stab_mask_2)

            # obtain ribbon labels here
            ribbonLabelString = ""

            uniques, indices = np.unique(
                okRoiAalLabels[commonNodes], return_inverse=True)
            counts = [np.sum(indices == k) for k in range(0, len(uniques))]

            ordering_big_to_small = np.argsort(counts)[::-1]

            counter = 0
            for index in ordering_big_to_small:
                unique = uniques[index]
                count = counts[index]
                if unique:
                    if count > ribbon_label_size_lim:
                        if counter < 3:
                            ribbonLabelString += unique + ", "
                            counter += 1
                        else:
                            ribbonLabelString += unique + "\n "
                            counter = 0
            ribbonLabelString = ribbonLabelString[:-2]
            ribbon_label_mat[i, j] = ribbonLabelString
            # return
    assert len(moduleSizes1) == len(mod_colors1)
    assert len(moduleSizes2) == len(mod_colors2)

    if mask_large_cluster_with_gray:
        mask_large_cluster_with_gray(mod_colors1, moduleSizes1)
        mask_large_cluster_with_gray(mod_colors2, moduleSizes2)

    fig = plt.figure(figsize=(10, 12))
    ax = fig.add_axes((0.05, 0.05, 0.90, 0.90))  # (111)

    alluvial.plot_alluvial(
        ax,
        ribbon_size_mat,
        ribbon_label_mat,
        mod_colors1,
        mod_colors2,
        ribbon_bglim=30,
        stable_ribbon_sizes_1=stableRibbonSizes1,
        stable_ribbon_sizes_2=stableRibbonSizes2,
        horizontal_pad_lr=0.0,
        vertical_pad_btw_modules=0.02,
        module_width=0.15
    )
    return fig


def mask_large_cluster_with_gray(mod_colors, mod_sizes):
    if len(mod_colors[0]) == 4:
        grey = (0.65, 0.65, 0.65, 1)
    else:
        grey = (0.65, 0.65, 0.65)
    mod_colors[mod_sizes > 0.5 * np.sum(mod_sizes)] = grey


def get_module_colors(sorted_clu_labels, n_cols=None):
    # indexing starts from zero
    if n_cols is None:
        n_cols = np.max(sorted_clu_labels) + 1
    cols = colorhelp.get_distinct_colors(n_cols)
    colors = cols[sorted_clu_labels]
    assert len(colors) == len(sorted_clu_labels)
    return colors


def add_colorbar(ax, vmin=0, vmax=1,
                 colormap="jet",
                 orientation="vertical",
                 cmapLabel="",
                 norm=None,
                 p_val_tick_vals=None,
                 p_val_tick_labels=None):
    """
    Add a colorbar to the axis.
    """
    if norm is None:
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb = mpl.colorbar.ColorbarBase(
        ax, cmap=colormap, norm=norm, orientation='vertical')
    cb.set_label(cmapLabel, rotation=90)
    if p_val_tick_vals:
        cb.set_ticks(p_val_tick_vals)
        cb.set_ticklabels(p_val_tick_labels)


def cluster_diff_plot_with_matrices(
        clustering, adj_mats, nodeColors, n1, paired_study, nIt=1e5,
        normcoords=None, suptitle="clusterplot", vmin=None, vmax=None,
        xlabel=None, invert_cmap=True,
        titles=None):
    # clustering = clustering[clustering!=settings.undef_clu_label]
    # clustering = communities.makeZeroIndexedClustering(clustering)

    movie_avg_adj_mat = np.average(adj_mats[:n1], axis=0)
    rest_avg_adj_mat = np.average(adj_mats[n1:], axis=0)
    abs_diff_mat = np.abs(movie_avg_adj_mat - rest_avg_adj_mat)

    t_val_diff_mat = measures.paired_t_value(
        np.array(adj_mats), len(adj_mats) / 2)

    fig = plt.figure(figsize=(12, 4))
    from matplotlib import gridspec
    gs = gridspec.GridSpec(1, 4, left=None, bottom=None, right=None, top=None,
                           wspace=None, hspace=None,
                           width_ratios=[1, 1, 1, 0.15])

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])
    ax3 = plt.subplot(gs[0, 2])

    if invert_cmap:
        cmap = cm.get_cmap("RdBu")
    else:
        cmap = cm.get_cmap('RdBu_r')


    # cols[len(cols) / 2] = (0.9, 0.9, 0.9, 1.0)

    # n_colors = 11  # pretty number a hack
    # tmp_colors = colormap(np.linspace(0, 1, n_colors))
    # colormap = mpl.colors.ListedColormap(tmp_colors)

    adj_mats = np.array(adj_mats)
    p_vals, mean_diffs = ptests.mean_difference_permtest(
        adj_mats, n1, len(adj_mats) - n1, paired_study, nIt)
    print p_vals

    p_vals_vec = p_vals[np.triu_indices_from(p_vals)]
    p_t, sig = corrections.fdr_bh(
        0.05, p_vals_vec, return_significant=True)


    if p_t < 0.001 or p_t > 0.005:
        p_val_limits = [0, 0.001, 0.05, 1, 1.95, 2 - 0.001, 2]
        # p_val_limits = [0, 0.001, 0.005, 0.01, 1.0, 1.99, 1.995, 1.999, 2]
        ticks = [str(i)
                 for i in [0, 0.001, 0.05, 1, 0.05, 0.001, 0]]
    else:
        p_val_limits = [0, 0.001, p_t, 0.05, 1, 1.95, 2 - p_t, 2 - 0.001, 2]
        # p_val_limits = [0, 0.001, 0.005, 0.01, 1.0, 1.99, 1.995, 1.999, 2]
        ticks = [str(i)
                 for i in [0, 0.001, str(p_t)[:6] + " FDR $<$ 0.05", 0.05, 1, 0.05, str(p_t)[:6] + " FDR $<$ 0.05", 0.001, 0]]
    assert len(p_val_limits) == len(ticks)
    # [0, 0.001, 0.005, 0.01, 1.0, 0.01, 0.005, 0.001, 0]]
    print p_vals_vec
    print np.sum(p_vals_vec <= p_t)
    print p_t
    print p_vals < p_t


    # debug###############
    assert (mean_diffs == movie_avg_adj_mat - rest_avg_adj_mat).all()
    one_sided_p_vals = p_vals / 2.0
    col_p_vals = 2 - ((1 + np.sign(mean_diffs)) - np.sign(mean_diffs) * p_vals)

    cols = cmap(np.array(np.linspace(0, 1, len(p_val_limits) - 1)))
    colormap, norm = colors.from_levels_and_colors(
        p_val_limits, cols, extend=u'neither')
    colorMatrix = colormap(norm(col_p_vals))

    max_ax12 = np.maximum(
        np.max(movie_avg_adj_mat), np.max(rest_avg_adj_mat)) * 1.2
    # weights = [1000, 3000, 10000, 30000, 50000]
    weights = [500, 1500, 5000, 15000, 50000]




    _hinton(ax1, movie_avg_adj_mat,
            clusterColors=nodeColors, xlabel=xlabel, max_weight=max_ax12, scale_weights=weights)
    _hinton(ax2, rest_avg_adj_mat,
            clusterColors=nodeColors, xlabel=xlabel, max_weight=max_ax12, scale_weights=weights)



    # xlim = ax1.get_xlim()
    # bbox = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    # width, height = bbox.width, bbox.height
    # max_esize = width / (xlim[1] - xlim[0])

    # ax_ascale = plt.subplot(gs[1, 0:2])
    # square_scale(ax_ascale, weights, max_width_in_inches=max_esize)
    # ax_dscale = plt.subplot(gs[1, 2])
    # exit()

    max_diff = np.max(abs_diff_mat)
    wscale_diff = [100, 300, 1000, 3000, 10000]
    _hinton(ax3, abs_diff_mat, clusterColors=nodeColors,
            colorMatrix=colorMatrix, xlabel=xlabel,
            scale_weights=wscale_diff,
            scale_box_color='white',
            max_weight=max_diff)  # max_diff)

    ax_cb = fig.add_axes([0.84, 0.3, 0.03, 0.47])

    add_colorbar(ax_cb, vmin=0, vmax=1, colormap=colormap,
                 orientation="vertical",
                 cmapLabel='p-value',
                 norm=norm,
                 p_val_tick_vals=p_val_limits,
                 p_val_tick_labels=ticks)

    if titles:
        if invert_cmap:
            ax1.set_title(titles[0], color=cols[-2])
            ax2.set_title(titles[1], color=cols[1])
        else:
            ax1.set_title(titles[0], color=cols[1])
            ax2.set_title(titles[1], color=cols[-2])
        ax3.set_title(titles[2])

        ax_cb.text(
            0.5, -0.05, r'\noindent \begin{center} More links \\in ' + titles[0].lower() + "\end{center}", va='top', ha='center', )
        ax_cb.text(
            0.5, 1.05, r'\noindent \begin{center} More links \\in ' + titles[1].lower() + "\end{center}", va='bottom', ha='center')

    fig.suptitle(suptitle)
    return fig


def corr_mat_diff_plot(
        corr_mats,
        mod_colors,
        n1, nIt,
        paired_study,
        normcoords=None,
        suptitle="clusterplot",
        vmin=None, vmax=None,
        xlabel=None, invert_cmap=True, titles=None):

    avg_adj_mat_1 = np.average(corr_mats[:n1], axis=0)
    avg_adj_mat_2 = np.average(corr_mats[n1:], axis=0)
    abs_diff_mat = np.abs(avg_adj_mat_1 - avg_adj_mat_2)

    fig = plt.figure(figsize=(12, 4))
    from matplotlib import gridspec
    gs = gridspec.GridSpec(1, 4, left=None, bottom=None, right=None, top=None,
                           wspace=None, hspace=None,
                           width_ratios=[1, 1, 1, 0.15])

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])
    ax3 = plt.subplot(gs[0, 2])

    if invert_cmap:
        cmap = cm.get_cmap("RdBu")
    else:
        cmap = cm.get_cmap('RdBu_r')

    # cols[len(cols) / 2] = (0.9, 0.9, 0.9, 1.0)

    # n_colors = 11  # pretty number a hack
    # tmp_colors = colormap(np.linspace(0, 1, n_colors))
    # colormap = mpl.colors.ListedColormap(tmp_colors)

    adj_mats = np.array(corr_mats)

    p_vals, mean_diffs = ptests.mean_difference_permtest(
        np.arctanh(adj_mats), n1, len(adj_mats) - n1, paired_study, nIt)

    p_vals_vec = p_vals[np.triu_indices_from(p_vals)]
    p_t, sig = corrections.fdr_bh(
        0.05, p_vals_vec, return_significant=True)

    p_val_limits = [0, 0.001, p_t, 0.05, 1, 1.95, 2 - p_t, 2 - 0.001, 2]
    # p_val_limits = [0, 0.001, 0.005, 0.01, 1.0, 1.99, 1.995, 1.999, 2]
    ticks = [str(i)
             for i in [0, 0.001, str(p_t)[:6] + " FDR $<$ 0.05", 0.05, 1, 0.05, str(p_t)[:6] + " FDR $<$ 0.05", 0.001, 0]]

    assert len(p_val_limits) == len(ticks)
    # [0, 0.001, 0.005, 0.01, 1.0, 0.01, 0.005, 0.001, 0]]
    print p_vals_vec
    print np.sum(p_vals_vec <= p_t)
    print p_t
    print p_vals < p_t

    # debug###############

    # assert (mean_diffs == avg_adj_mat_1 - avg_adj_mat_2).all()

    # one_sided_p_vals = p_vals / 2.0
    col_p_vals = 2 - ((1 + np.sign(mean_diffs)) - np.sign(mean_diffs) * p_vals)

    cols = cmap(np.array(np.linspace(0, 1, len(p_val_limits) - 1)))
    colormap, norm = colors.from_levels_and_colors(
        p_val_limits, cols, extend=u'neither')
    color_mat = colormap(norm(col_p_vals))

    max_ax12 = np.maximum(
        np.max(avg_adj_mat_1), np.max(avg_adj_mat_2)) * 1.2
    # weights = [1000, 3000, 10000, 30000, 50000]

    corr_cmap = cm.get_cmap('RdBu_r')

    max_avg_corr = np.maximum(np.max(avg_adj_mat_1), np.max(avg_adj_mat_2))
    _hinton(ax1, np.ones(avg_adj_mat_1.shape),
            clusterColors=mod_colors, xlabel=xlabel,
            colorMatrix=corr_cmap((avg_adj_mat_1)/max_avg_corr),
            max_weight=max_ax12
    )

    _hinton(ax2, np.ones(avg_adj_mat_2.shape),
            clusterColors=mod_colors, xlabel=xlabel,
            colorMatrix=corr_cmap((avg_adj_mat_2)/max_avg_corr),
            max_weight=max_ax12
    )

    if invert_cmap:
        ax1.set_title(titles[0], color=cols[-2])
        ax2.set_title(titles[1], color=cols[1])
    else:
        ax1.set_title(titles[0], color=cols[1])
        ax2.set_title(titles[1], color=cols[-2])

    _hinton(ax3, np.ones(abs_diff_mat.shape),
            clusterColors=mod_colors,
            colorMatrix=color_mat,
            xlabel=xlabel,
            scale_box_color='white')
    ax3.set_title(titles[2])

    ax_cb = fig.add_axes([0.84, 0.3, 0.03, 0.47])

    add_colorbar(ax_cb, vmin=0, vmax=1, colormap=colormap,
                 orientation="vertical",
                 cmapLabel='p-value',
                 norm=norm,
                 p_val_tick_vals=p_val_limits,
                 p_val_tick_labels=ticks)
    if titles:
        ax_cb.text(
            0.5, -0.05, r'\noindent \begin{center} More correlated \\in ' + titles[0].lower() + "\end{center}", va='top', ha='center', )
        ax_cb.text(
            0.5, 1.05, r'\noindent \begin{center} More correlated \\in ' + titles[1].lower() + "\end{center}", va='bottom', ha='center')

    fig.suptitle(suptitle)
    return fig


def _hinton(ax, sizeMatrix, clusterColors=None, colorMatrix=None,
            xlabel=True, max_weight=None, bg_color=None, box_color=None,
            scale_weights=None, scale_box_color=None):
    """Draw Hinton diagram for visualizing a weight matrix.

    Parameters
    ----------
    ax : matplotlib axes to plot the hinton diagram to
    sizeMatrix : 2d numpy matrix of shape (n,n)
        describes the areas of the squares in the hinton diagram
    clusterColors : list-like
        list of colors for each cluster
    colorMatrix : 2d numpy array of shape (n,n)
        each element should correspond to a matplotlib color
    xlabel : str
        what to use for labeling for x and y axis
    max_weight : float
        the maximum "size of a block"
    bg_color: matplotlib color
        'gray'
    box_color: float between 0.0 and 1.0
        shades of grey at the moment

    Returns
    -------
    None
    """
    ax = ax if ax is not None else plt.gca()

    if bg_color is None:
        bg_color = 'white'
    if box_color is None:
        box_color = 0.0

    if max_weight is None:
        max_weight = np.max(sizeMatrix) * 1.1
    if colorMatrix is None:
        x, y = sizeMatrix.shape
        colorMatrix = np.ones((x, y, 3)) * box_color

    sizeMatrix = sizeMatrix[::-1, ::-1]
    colorMatrix = colorMatrix[::-1, ::-1]
    if clusterColors is not None:
        clusterColors = clusterColors[::-1]

    ax.patch.set_facecolor(bg_color)
    ax.set_aspect('equal', 'box')

    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())

    for (i, j), w in np.ndenumerate(sizeMatrix):
        color = colorMatrix[i, j]
        side_length = np.sqrt(w / float(max_weight))
        # print i, j, w, side_length
        rect = plt.Rectangle(
            [i - side_length / 2, j - side_length /
                2], side_length, side_length,
            facecolor=color, edgecolor=color)
        ax.add_patch(rect)

    if clusterColors is not None:
        for i, color in enumerate(clusterColors):
            rect = plt.Rectangle(
                [i - 1. / 2., - 3 / 2.], 1, 1, facecolor=color,
                edgecolor=color)
            ax.add_patch(rect)
            rect = plt.Rectangle(
                [- 3. / 2., i - 1 / 2.], 1, 1, facecolor=color,
                edgecolor=color)
            ax.add_patch(rect)

        ax.add_patch(
            plt.Rectangle(
                (-0.65, -1.5),
                0.15,
                len(clusterColors) + 1,
                facecolor="k",
                edgecolor="k",
                lw=1
            )
        )
        ax.add_patch(
            plt.Rectangle(
                (-1.5, -0.65),
                len(clusterColors) + 1,
                0.15,
                facecolor="k",
                edgecolor="k",
                lw=1
            )
        )

        x, y = sizeMatrix.shape
        ax.set_xlim([-1.5, x - 0.5])
        ax.set_ylim([-1.5, y - 0.5])
        ax.set_xlabel(xlabel)
        ax.xaxis.set_label_position('top')
        ax.set_ylabel(xlabel)
    else:
        x, y = sizeMatrix.shape
        ax.set_xlim([-0.5, x - 0.5])
        ax.set_ylim([-0.5, y - 0.5])

    # scale
    if scale_weights is not None:
        ax.set_axis_off()
        # around rectangle
        ax.add_patch(plt.Rectangle(
            [- 3. / 2., - 3 / 2.],
            len(clusterColors) + 1,
            len(clusterColors) + 1,
            facecolor=[0, 0, 0, 0],
            edgecolor="k",
            lw=1)
        )
        ax.set_xlim(-4.5, len(clusterColors) + 1.1)
        ax.set_ylim(-1.6, len(clusterColors) + 1.1)

        y_locs = np.linspace(len(clusterColors) - 2, 2, len(scale_weights))
        x_center = -3

        if scale_box_color is None:
            scale_box_color = 'k'

        for i, w in enumerate(scale_weights):
            side_length = np.sqrt(w / float(max_weight))
            rect = plt.Rectangle(
                [x_center - side_length / 2, y_locs[i] - side_length /
                 2], side_length, side_length,
                facecolor=scale_box_color, edgecolor='k')
            ax.add_patch(rect)
            ax.text(x_center - 0.9, y_locs[i], str(w),
                    rotation=00, ha='right', va='center')
        ax.text(
            x_center - 1.5, y_locs[-1] - 1.7,
            r"\noindent Number \\ \noindent of links",
            ha='center', va='center')
    else:
        ax.autoscale_view()
    ax.invert_yaxis()


def get_hash_name(string_list):
    """
    Hashs a list of strings to create a semi-temporary filename
    """
    string = ""
    for s in string_list:
        string += str(s)
    return hashlib.sha1(string).hexdigest()


class RecomputeException(Exception):
    pass


def get_clu_net_adj_mats(cfg,
                         corr_mat_fnames,
                         clustering,
                         recompute=True,
                         weighted=False):
    clunet_adj_mat_fname = (cfg['outdata_dir'] + "/" +
                            get_hash_name([os.getcwd()] +
                                          corr_mat_fnames +
                                          list(clustering) +
                                          [cfg['density']]) +
                            ".pkl")

    try:
        if recompute:
            print "recompute!"
            raise RecomputeException("Recompute!")
        # check if already computed
        out_dict = dataio.load_pickle(clunet_adj_mat_fname)
    except (RecomputeException, IOError) as e:
        clu_net_adj_mats = []
        for fname in corr_mat_fnames:
            print fname
            # for presenting the adj. matrix
            cont_zero_indexed_clu = communities.make_zero_indexed_clustering(
                clustering)
            ok_nodes = dataio.get_ok_nodes(cfg['blacklist_fname'])
            net = netgen.get_graph_from_bare_data(
                fname, cfg['blacklist_fname'], cfg['density'],
                include_mst=cfg['include_mst'], weighted=False
            )
            clu_net_adj_mat = communities.comp_cluster_adj_mat(
                cont_zero_indexed_clu[ok_nodes],
                net
            )
            clu_net_adj_mats.append(clu_net_adj_mat)

        out_dict = {"clu_net_adj_mats": clu_net_adj_mats}
        print clunet_adj_mat_fname
        dataio.save_pickle(clunet_adj_mat_fname, out_dict)
    return out_dict["clu_net_adj_mats"]


def comp_and_viz_cluster_diff_matrices_from_file(
        fname,
        cfg,
        corr_mat_fnames,
        n1,
        **kwargs):
    """
    Simple wrapper for :py:func:`comp_and_viz_cluster_diff_matrices`
    to use a filename instead of a network partition.

    Parameters
    ----------
    fname : str
        path to file containing the clustering

    For others parameters and kwargs, see
    :py:func:`comp_and_viz_cluster_diff_matrices`
    """

    return comp_and_viz_cluster_diff_matrices(
        dataio.load(fname)[settings.louvain_cluster_tag],
        cfg,
        corr_mat_fnames,
        n1,
        **kwargs
    )


def comp_and_viz_cluster_diff_matrices(clustering,
                                       cfg,
                                       corr_mat_fnames,
                                       n1,
                                       vmin=-5,
                                       vmax=5,
                                       suptitle="",
                                       recompute=True,
                                       mod_colors=None,
                                       invert_cmap=False,
                                       return_cg_mats=False,
                                       titles=None):
    """
    Parameters
    ----------
    clustering : list, numpy array
        a cluster assignment membership list, with blacklisted nodes
    cfg : dict
        brainnets config dict
    corr_mat_fnames : list of str
        the corr_mat_fnames
    n1 : int
        number of elements in the first group
    vmin : float
        vmin for the colormap
    vmax : float
        vmax for the colormap
    suptitle : str
        title for the str
    recompute : bool
        recompute the correlation matrix

    Returns
    -------
    fig : matplotlib figure
        the figure with

    """

    clustering = np.array([int(label) for label in clustering])
    clu_labels = np.sort(
        np.setdiff1d(np.unique(clustering), [settings.undef_clu_label]))

    clu_net_adj_mats = get_clu_net_adj_mats(
        cfg,
        corr_mat_fnames,
        clustering,
        recompute=recompute,
        weighted=False)

    print clu_net_adj_mats[0].shape

    if mod_colors is None:
        mod_colors = get_module_colors(clu_labels)
    mod_sizes = communities.get_module_sizes(
        clustering[clustering != settings.undef_clu_label])

    mask_large_cluster_with_gray(mod_colors, mod_sizes)

    fig = cluster_diff_plot_with_matrices(
        clustering,
        clu_net_adj_mats,
        mod_colors,
        n1,
        cfg['paired'],
        "all",
        suptitle=suptitle,
        invert_cmap=invert_cmap,
        titles=titles
    )
    if return_cg_mats:
        return fig, clu_net_adj_mats
    else:
        return fig


def plotAllCommunityVisualizations(
        clu_1_fname, clu_2_fname, cfg,
        corr_mat_fnames=[], density=None, n1=None,
        nIt=None, tValVmin=-6, tValVmax=6, ribbon_label_size_lim=15,
        recompute=False, stab_mask_1=None, stab_mask_2=None,
        plot_alluvial=True, plot_on_brain=True, plot_diff_nets=True,
        plot_diff_matrices=True):
    """
    Plots different community structure visualizations.
    """
    print "plot all community visualizations needs updating"
    return None
    figsToReturn = []

    # plot alluvial
    if plot_alluvial:
        aFig = plot_alluvial_diagram(
            clu_1_fname, clu_2_fname, cfg[
                'node_info_fname'], cfg['blacklist_fname'],
            ribbon_label_size_lim, stab_mask_1=stab_mask_1,
            stab_mask_2=stab_mask_2)
        figsToReturn.append(aFig)

    # plot the actual communities
    clu1 = dataio.loadPickle(clu_1_fname)[settings.louvain_cluster_tag]
    clu2 = dataio.loadPickle(clu_2_fname)[settings.louvain_cluster_tag]
    clu1, clu2, n_clu_1, n_clu_2 = _get_matched_renumbered_communities_for_viz(
        clu1, clu2)
    # get colors here?

    if plot_on_brain:
        brain_fig1 = viz_com_structure_using_slices(
            clu1, cfg['node_info_fname'], cfg['blacklist_fname'])
        figsToReturn.append(brain_fig1)
        brain_fig2 = viz_com_structure_using_slices(
            clu2, cfg['node_info_fname'], cfg['blacklist_fname'])
        figsToReturn.append(brain_fig2)

    # if plot_diff_nets:
    # all the tricky [0,1] stuff due to oords=None?
    #     if n_clu_1 > n_clu_2:
    #         clusInOrder = [0, 1]
    #     else:
    #         clusInOrder = [1, 0]

    #     clus = [clu1, clu2]
    #     cluIndexToFig = {}
    #     for cluIndex in clusInOrder:
    #         clu = clus[cluIndex]
    #         f, normcoords = computeAndVizClusterDiffNetwork(
    #             clu, corr_mat_fnames, cfg['blacklist_fname'], density, n1,
    #             paired_study, nIt,
    #             vmin=-5, vmax=5, suptitle="", normcoords=None,
    #             save_pickleData=True, recompute=recompute)
    #         cluIndexToFig[cluIndex] = f

    #     for clu in [0, 1]:
    #         figsToReturn.append(cluIndexToFig[clu])

    if plot_diff_matrices:
        if recompute and plot_diff_nets:
            recompute = False
        for clu in [clu1, clu2]:
            figsToReturn.append(
                comp_and_viz_cluster_diff_matrices(
                    clu, corr_mat_fnames, cfg['blacklist_fname'], density, n1,
                    cfg['paired'], nIt,
                    vmin=-5, vmax=5, suptitle="",
                    save_pickleData=True, recompute=recompute))
    return figsToReturn


def _convert_igraph_2_netpython(igraph_net, add_weight_to_zero_links=True):
    """
    Creates a pynet.SymmNet instance from a given adjacency matrix.
    Nodelabels are obtained from the 'names' parameter.

    **Use this function only for layouting purposes!**
    """
    from netpython import pynet
#   assert len(igraph_net.vs) == len(labels)
    net = pynet.SymmNet(len(igraph_net.vs))

    if igraph_net.summary()[9] == "W":
        weights = np.array(igraph_net.es["weight"])
    else:
        weights = np.ones(len(igraph_net.es))

    for i, edge in enumerate(igraph_net.es):
        #        if weights[i] == 0:
        if add_weight_to_zero_links:
            net[edge.source, edge.target] = (weights[i] +
                                             0.01 * np.max(weights))  # used only for layout
        else:
            net[edge.source, edge.target] = weights[i]
    return net

# OLD STUFF THAT WAS USED FOR PLOTTING DIFF MAT AS A NETWORK:
#
# def getNewNormCoords(oldClu, oldNormCoords, newClu):
#     oldCluLabels = np.sort(np.setdiff1d(oldClu, [settings.undef_clu_label]))
#     newCluLabels = np.sort(np.setdiff1d(newClu, [settings.undef_clu_label]))
#     oldNormCoordsCounter = 0
#     newNormCoords = []
#     for currentLabel in range(0, np.max(newClu) + 1):
#         if currentLabel in oldCluLabels:
#             if currentLabel in newCluLabels:
# print currentLabel, "cl", oldNormCoordsCounter, "oncc"
#                 newNormCoords.append(oldNormCoords[oldNormCoordsCounter])
#             oldNormCoordsCounter += 1
#         else:
#             if currentLabel in newCluLabels:
#                 newNormCoords.append(0.1 + 0.8 * np.random.random(2))
#     return newNormCoords

# def computeAndVizClusterDiffNetwork(clustering, corr_mat_fnames,
#                                     blacklist_fname, density, n1,
#                                     pairedStudy,
#                                     nIt, vmin=-5, vmax=5, suptitle="",
#                                     normcoords=None, save_pickleData=True,
#                                     recompute=True):
#     """Computes and plots the cluster network plot for a given clustering and
#     correlation matrices.

#     Parameters
#     ----------
#     clustering : [-1,2,3,4,1,2,....,1,1,1,2,2,1,-1,1]
#     corr_mat_fnames : names of the correlation matrices
#     density : network desnity to be used in corr.mat -> net
#     n1 : number of subs in the first group
#     pairedStudy : True/False
#     nIt : e.g. 10 or "all" (if pairedStudy==True)
#     vmin, vmax: limits for the colorbar
#     suptitle : e.g. "movie clus plaa"
#     normcoords : coordinates for the
#     """

#     cluLabels = np.sort(
#         np.setdiff1d(np.unique(clustering), [settings.undef_clu_label]))

#     clu_net_adj_mats = get_clu_net_adj_mats(
#         corr_mat_fnames, clustering, density, blacklist_fname,
#         recompute=recompute, weighted=False, includeMST=True
#     )
#     nodeColors = get_module_colors(cluLabels)
#     if normcoords is not None:
#         newcoords = []
#         for i in np.sort(np.unique(clustering)):
#             newcoords.append(normcoords[i])
#         normcoords = np.array(newcoords)
#     fig, normcoords = clusterDiffPlot(
#         clustering, clu_net_adj_mats, nodeColors, n1,
#         pairedStudy, "all", normcoords=normcoords,
#         suptitle=suptitle, vmin=-5, vmax=5)
#     return fig, normcoords

# def clusterDiffPlot(
#         clustering, adjMatrices, nodeColors, n1, pairedStudy, nIt=1e5,
#         normcoords=None, suptitle="clusterplot", vmin=None, vmax=None):

#     def weight2ew(weights, min_w=None, max_w=None):
#         if min_w is None:
#             min_w = np.min(weights)
#         if max_w is None:
#             max_w = np.max(weights)
#         edge_widths = 0.01 + np.array(weights) / (max_w) * 10
#         return edge_widths

#     clustering = clustering[clustering != settings.undef_clu_label]
#     clustering = communities.makeZeroIndexedClustering(clustering)

#     cluSizes = np.zeros(len(nodeColors))
#     for i, clu in enumerate(np.sort(np.unique(clustering))):
#         cluSizes[i] = np.sum(clustering == clu)

#     cluNodeSizes = np.sqrt(cluSizes) / np.sqrt(np.max(cluSizes)) * 20

#     cluLabels = ["" for i in range(len(nodeColors))]

# compute circular coords?

#     fig = plt.figure(figsize=(24, 8))
#     ax1 = fig.add_axes([0.00, 0.0, 0.28, 1])
#     ax2 = fig.add_axes([0.30, 0.0, 0.28, 1])
#     ax3 = fig.add_axes([0.60, 0.0, 0.28, 1])
#     ax_cb = fig.add_axes([0.94, 0.10, 0.02, 0.8])

# avgCluNetAdjMatrixMode1 = np.average(adjMatrices[:nMode1], axis=0)
#     movieAverageAdjMatrix = np.average(adjMatrices[:n1], axis=0)
#     restAverageAdjMatrix = np.average(adjMatrices[n1:], axis=0)

#     graph1, weights1 = netgen.makeFullWeightedNetFromAdjMatrix(
#         movieAverageAdjMatrix, return_weights=True)
#     graph2, weights2 = netgen.makeFullWeightedNetFromAdjMatrix(
#         restAverageAdjMatrix, return_weights=True)
#     minw = np.min([weights1, weights2])
#     maxw = np.max([weights1, weights2])

#     edge_widths = weight2ew(weights1, min_w=minw, max_w=maxw)
#     _, normcoords = plot_network(
#         igraph_net=graph1, nodeColors=nodeColors, axes=ax1,
#         normcoords=normcoords,
#         nodeLabels=cluLabels, nodeSizes=cluNodeSizes,
#         edgeWidths=edge_widths)

#     edge_widths = weight2ew(weights2, min_w=minw, max_w=maxw)

#     _, normcoords = plot_network(
#         igraph_net=graph2, nodeColors=nodeColors, axes=ax2,
#         normcoords=normcoords,
#         nodeLabels=cluLabels, nodeSizes=cluNodeSizes,
#         edgeWidths=edge_widths)

# make the diff plot
#     diffCluNetAdjMatrix = movieAverageAdjMatrix - restAverageAdjMatrix

# compute stats for edges
# allWeights = []  # number of links between modulse
# allNCluInternalLinks = []  # number of links withing module
#     for cluNetAdjMatrix in adjMatrices:
#         allWeights.append(
#             cluNetAdjMatrix[np.triu_indices_from(cluNetAdjMatrix, 1)])
#         allNCluInternalLinks.append(np.diag(cluNetAdjMatrix))
# compute stats for nodes

#     tValues = measures.paired_t_value(
#         np.array(allWeights), len(allWeights) / 2)
#     tValues[np.isnan(tValues)] = 0
#     tValues_nodes = measures.paired_t_value(
#         np.array(allNCluInternalLinks), len(allNCluInternalLinks) / 2)
#     tValues_nodes[np.isnan(tValues_nodes)] = 0
#     pValues, meanDiffs = ptests.mean_difference_permtest(
#         np.array(allWeights), n1, n1, True, "all")
#     pValues_node, meanDiffs = ptests.mean_difference_permtest(
#         np.array(allNCluInternalLinks), n1, n1, True, "all")

#     linkLabels = []
#     for pValue in pValues:
#         if pValue < 0.05:
#             linkLabels.append(str(np.log10(pValue))[:4])
#         else:
#             linkLabels.append("")

#     nodeLabels = []
#     for pValue in pValues_node:
#         if pValue < 0.05:
#             nodeLabels.append(str(np.log10(pValue))[:4])
#         else:
#             nodeLabels.append("")

# linkLabels = [str(np.log10(pValue))[:4] for pValue in pValues]
# cluNodeLabels = [str(np.log10(pValue))[:4] for pValue in pValues_node]

#     absmaxTVal = np.max(
#         [np.max(np.abs(tValues)), np.max(np.abs(tValues_nodes))])
#     colormap = cm.get_cmap("RdBu_r")
#     linkColors = colormap(0.5 * (tValues / absmaxTVal + 1))
#     nodeTValColors = colormap(0.5 * (tValues_nodes / absmaxTVal + 1))

#     graph, weights = netgen.makeFullWeightedNetFromAdjMatrix(
#         diffCluNetAdjMatrix, return_weights=True)
#     edge_widths = weight2ew(np.abs(weights), min_w=minw, max_w=maxw)
#     plot_network(
#         igraph_net=graph, nodeColors=nodeTValColors, normcoords=normcoords,
#         axes=ax3, edgeWidths=edge_widths,
#         edgeColors=linkColors, edgeLabels=linkLabels, edgeLabelSize=10,
#         nodeSizes=cluNodeSizes, nodeLabels=nodeLabels)

#     if vmin is None:
#         vmin = -1 * absmaxTVal
#     if vmax is None:
#         vmax = absmaxTVal
#     add_colorbar(ax_cb, vmin=vmin, vmax=vmax, colormap=colormap,
#                  orientation="vertical",
#                  cmapLabel=settings.get_prop_tex_name(settings.tval_tag))

#     fig.suptitle(suptitle)
#     return fig, normcoords


# def plot_network(netpython_net=None,
#                  igraph_net=None,
#                  axes=None,
#                  normcoords=None,
#                  nodeLabels=None,
#                  edgeLabels=None,
#                  edgeLabelSize=None,
#                  nodeColors=None,
#                  nodeSizes=None,
#                  nodeLabelSizes=None,
#                  nodeLabelOffsets=None,
#                  edgeWidths=None,
#                  edgeColors=None,
#                  defaultNodeSize=10,
#                  ax_offset=0.1):
#     """
#     Plot a graph (igraph or netpython)
#     Give either a netpython_net or an igraph_net
#     (an igraph_net will be converted to a netpython_net internally)

#     Parameters
#     ----------
#     netpython_net : netpython network
#         to be plotted
#     igraph_net : igraph.Graph
#         instance to be plotted, used if netpython_net not specified
#     axes : matplotlib axes object
#         the axes to plot to
#     normcoords : dict
#         normalized coordinates of the network which can be used for plotting
#     nodeLabels : list
#         labels of the nodes (if igraph network)
#     nodeColors : list
#         colors of the nodes
#     nodeSizes : list
#         specify every node's size by giving a list of nodes
#     edgeWidths : list, optional
#         Defaults to 2
#     edgeColors : list
#         Defaults to gray ("0.5")
#     defaultNodeSize: float
#         the default node size if nodeSizes is not defined

#     Returns
#     -------:
#     fig : a matplotlib figure
#     normcoords :
#         the normalized coordinates of the nodes
#     """
# netpython modules needed
#     try:
#         from netpython import visuals
#     except:
#         print "you need the netpython lib developed"
#         print "at BECS networks group for this"
#         print "try to get it from somewhere... "
#         return

#     if netpython_net is None:
#         if igraph_net is None:
#             print "You should give at least one network to plot"
#         else:
# assume the igraph_net exist -> create a corresponding netpython
# network to be able to use functions from
#             netpython_net = __convert_igraph_2_netpython(igraph_net)

#     if normcoords is None:
#         coords = visuals.Himmeli(netpython_net).getCoordinates()
#         coordvals = np.array(coords.values())[:, :2]
#         [xmax, ymax] = np.max(coordvals, 0)
#         [xmin, ymin] = np.min(coordvals, 0)
#         xran = xmax - xmin
#         yran = ymax - ymin
#         normcoords = {}
#         for key in coords:
#             normx = coords[key][0] / xran * (1 - 2 * ax_offset) + ax_offset
#             normy = coords[key][1] / yran * (1 - 2 * ax_offset) + ax_offset
#             normcoords[key] = [normx, normy]

#     if nodeSizes is None:
#         nodeSizes = [defaultNodeSize for i in range(len(normcoords))]

#     fig = None
#     if axes is None:
#         fig = plt.figure(figsize=(6, 6))
#         axes = fig.add_subplot(111)

# plot edges
#     if edgeWidths is None:
#         edgeWidths = [2 for i in range(len(netpython_net.edges))]
#     if edgeColors is None:
#         edgeColors = ["0.5" for i in range(len(netpython_net.edges))]

# to calibrate edge_widths etc..
# edgelist = list(netpython_net.edges)  # , dtype=np.int32)
#     nEdges = len(edgelist)
#     edgelist = sorted(
#         edgelist, key=lambda column_entry:
#         column_entry[0] * nEdges + column_entry[1])

#     for c, [i, j, w] in enumerate(edgelist):
#         xcoords = [normcoords[i][0], normcoords[j][0]]
#         ycoords = [normcoords[i][1], normcoords[j][1]]
#         axes.plot(
#             xcoords, ycoords, color=edgeColors[c], ls='-', lw=edgeWidths[c])
#         if edgeLabels is not None:
#             if edgeLabelSize:
#                 size = edgeLabelSize
#             else:
#                 size = max(edgeWidths[c], defaultNodeSize * 2)
#             axes.annotate(edgeLabels[c], (np.mean(xcoords),
#                           np.mean(ycoords)),
#                           textcoords='offset points',
#                           xytext=(0, 0),
#                           color="k",
#                           va = "center",
#                           ha = "center",
#                           size=size)

# plot nodes and nodelabels
#     for nodeIndex in netpython_net:
# print nodeIndex, len(nodeSizes)
#         nodeSize = nodeSizes[nodeIndex]
#         color = "b"
#         if nodeColors is not None:
#             color = nodeColors[nodeIndex]
#         x, y = normcoords[nodeIndex][:2]
#         axes.plot([x], [y], "o", color=color, markersize=nodeSize)
#         if nodeLabels is not None:
#             label = nodeLabels[nodeIndex]
#         else:
#             label = nodeIndex
#         if nodeLabelSizes:
#             nodeLabelSize = nodeLabelSizes[nodeIndex]
#         else:
#             nodeLabelSize = max(nodeSize, defaultNodeSize * 2)

#         if nodeLabelOffsets:
#             nodeLabelOffset = nodeLabelOffsets[nodeIndex]
#         else:
#             nodeLabelOffset = (nodeSize / 2, nodeSize / 2)

#         axes.annotate(
#             label, (normcoords[nodeIndex][0], normcoords[nodeIndex][1]),
#             textcoords='offset points',
#             xytext=nodeLabelOffset,
#             va='center',
#             ha='center',
#             color="k",
#             size=nodeLabelSize)
# max(nodeSize, defaultNodeSize*2))

#     axes.set_frame_on(False)
#     axes.set_xticklabels([])
#     axes.xaxis.set_ticks_position('none')
#     axes.set_yticklabels([])
#     axes.yaxis.set_ticks_position('none')
#     axes.set_xlim([0, 1])
#     axes.set_ylim([0, 1])
#     return fig, normcoords
