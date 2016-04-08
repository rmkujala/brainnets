# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 10:08:05 2013

@author: rmkujala
"""
# third party
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib import rc
# brainnets imports
from brainnets import dataio, gencomps, settings, aux, config
from brainnets import genplots
from brainnets import fname_conventions as fnc

from verkko.permtests import bootstrap
from verkko.permtests import ptests
from verkko.permtests import measures

# matplotlib rc parameters
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern serif']})
rc('legend', fontsize=12)


def plot_pooled_correlation_dists_by_condition(cfg):
    """
    Plot the pooled (ie. combined, pooled) correlation distributions
    separately for each of the conditions (each list in fname_group_list)
    corresponds to a condition.

    Parameters
    ----------
    cfg : dict
        the user-specified config dictionary.
        The following keys are required::

            "group_1_mat_fnames"
            "group_2_mat_fnames"
            "group_1_color"
            "group_2_color"
            "group_1_label"
            "group_2_label"

    Returns
    -------
    fig : the matplotlib figure object
    """
    # CFGCHANGE?
    config.require(cfg, ["group_1_mat_fnames", "group_2_mat_fnames",
                         "group_1_color", "group_2_color", "group_1_label",
                         "group_2_label"])
    fname_group_list = [cfg["group_1_mat_fnames"], cfg["group_2_mat_fnames"]]
    colors = [cfg["group_1_color"], cfg["group_2_color"]]
    labels = [cfg["group_1_label"], cfg["group_2_label"]]

    n_bins = 100
    corr_bins, corr_bin_centers = aux.get_lin_bins(n_bins, -1, 1.)
    bin_counts = [np.zeros(n_bins) for f in fname_group_list]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, fname_group in enumerate(fname_group_list):
        for fname in fname_group:
            flat_corr_mat = \
                dataio.get_blacklist_filtered_and_flattened_adj_mat(
                    fname, cfg["blacklist_fname"]
                )
            bin_counts[i] += aux.get_bin_counts(flat_corr_mat, corr_bins)
        # normalize
        bin_counts[i] = bin_counts[i] / \
            (np.sum(bin_counts[i] * (corr_bins[1] - corr_bins[0])))
        ax.plot(corr_bin_centers, bin_counts[
                i], color=colors[i], label=labels[i])

    ax.set_xlabel(settings.get_prop_tex_name(settings.correlation_tag))
    ax.set_ylabel(r"Probability density P(c)")
    ax.legend(loc=0)
    fig.savefig(cfg['outdata_dir'] + "pooledCorrDists.pdf",
                format="pdf", bbox_inches='tight')
    return fig


def plot_corr_dists_subjectwise(cfg):
    """
    Plotting the individual correlation profiles to see the variations.

    Parameters
    ----------
    cfg : dict
        the user-specified config dictionary.
        The following keys are required::

            "group_1_mat_fnames"
            "group_2_mat_fnames"
            "group_1_color"
            "group_2_color"
            "group_1_label"
            "group_2_label"
            "outdata_dir"

    Returns
    -------
    fig : the matplotlib figure object
    """
    fname_group_list = [cfg["group_1_mat_fnames"], cfg["group_2_mat_fnames"]]
    colors = [cfg["group_1_color"], cfg["group_2_color"]]
    labels = [cfg["group_1_label"], cfg["group_2_label"]]
    rc('legend', fontsize=8)
    fig = plt.figure()
    n_bins = 100
    n = len(fname_group_list[0])
    ax_x, ax_y = _get_subplot_x_y(n)
    for i, _ in enumerate(fname_group_list[0]):  # now n1 is n2)
        ax = fig.add_subplot(ax_y, ax_x, i + 1)
        for j, fname_group in enumerate(fname_group_list):
            corr_bins, corr_bin_centers = aux.get_lin_bins(n_bins, -1, 1.)
            flat_corr_mat = \
                dataio.get_blacklist_filtered_and_flattened_adj_mat(
                    fname_group[i], cfg["blacklist_fname"]
                )
            bin_counts = aux.get_bin_counts(flat_corr_mat, corr_bins)
            # normalize
            bin_counts = bin_counts / \
                (np.sum(bin_counts * (corr_bins[1] - corr_bins[0])))
            ax.plot(corr_bin_centers, bin_counts, color=colors[
                    j], label=labels[j] + "\_" + str(i))
#        ax.set_xlabel(settings.get_prop_tex_name("corr"))
#        ax.set_ylabel(r"Probability density P(c)")
        ax.legend(loc=0)
    plt.tight_layout()
    fig.savefig(
        cfg["outdata_dir"] + "individualCorrDistsSubjectwise.pdf",
        format="pdf", bbox_inches='tight')
    return fig


def plot_individual_correlation_dists(cfg):
    """
    Plotting the individual correlation profiles to see the variations.

    Parameters
    ----------
    cfg : dict
        the user-specified config dictionary.
        The following keys are required::

            "group_1_mat_fnames"
            "group_2_mat_fnames"
            "group_1_color"
            "group_2_color"
            "group_1_label"
            "group_2_label"
            "outdata_dir"

    Returns
    -------
    fig : the matplotlib figure object
    """
    # CFGCHANGE?
    config.require(cfg, )
    fname_group_list = [cfg["group_1_mat_fnames"], cfg["group_2_mat_fnames"]]
    colors = [cfg["group_1_color"], cfg["group_2_color"]]
    labels = [cfg["group_1_label"], cfg["group_2_label"]]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    n_bins = 100
    corr_bins, corr_bin_centers = aux.get_lin_bins(n_bins, -1, 1.)

    for i, fname_group in enumerate(fname_group_list):
        for fname in fname_group:
            flat_corr_mat = \
                dataio.get_blacklist_filtered_and_flattened_adj_mat(
                    fname, cfg['blacklist_fname']
                )
            bin_counts = aux.get_bin_counts(flat_corr_mat, corr_bins)
        # normalize
            bin_counts = bin_counts / \
                (np.sum(bin_counts * (corr_bins[1] - corr_bins[0])))
            ax.plot(corr_bin_centers, bin_counts, color=colors[i],
                    label=labels[i])

    ax.set_xlabel(settings.get_prop_tex_name("corr"))
    ax.set_ylabel(r"Probability density P(c)")
    ax.legend(loc=0)
    fig.savefig(cfg['outdata_dir'] +
                "individualCorrDists.pdf", format="pdf", bbox_inches='tight')
    return fig


def plot_link_dist_probs_by_condition(cfg):
    """
    Plots the link distance PDFs pooled by condition.

    Parameters
    ----------
    cfg : dict
        the user-specified config dictionary.
        The following keys are required::

            "group_1_mat_fnames"
            "group_2_mat_fnames"
            "group_1_color"
            "group_2_color"
            "group_1_label"
            "group_2_label"
            "outdata_dir"

    Returns
    -------
    None
    """
    # CFGCHANGE?
    config.require(
        cfg, ["group_1_mat_fnames",
              "group_2_mat_fnames",
              "group_1_color",
              "group_2_color",
              "group_1_label",
              "group_2_label",
              "outdata_dir",
              "density"]
    )
    fname_group_list = [cfg["group_1_mat_fnames"], cfg["group_2_mat_fnames"]]
    colors = [cfg["group_1_color"], cfg["group_2_color"]]
    labels = [cfg["group_1_label"], cfg["group_2_label"]]

    # ns = []
    data_example = dataio.load_pickle(
        fnc.get_ind_fname(fname_group_list[0][0],
                          cfg, settings.link_distance_tag))
    densities = data_example[settings.config_tag][settings.densities_tag]

    # for individual plots:
    # figs = [plt.figure() for p in densities]
    # axs = [fig.add_subplot(111) for fig in figs]

    group_distances = []
    for i, fname in enumerate(fname_group_list):
        # ns.append(len(fname))
        distances = [np.array([]) for p in densities]
        for fname in fname:
            print fname
            data = dataio.load_pickle(
                fnc.get_ind_fname(fname, cfg, settings.link_distance_tag))
            p_f_dists = data[settings.link_distance_tag]
            for j, p in enumerate(densities):
                # print p
                distances[j] = np.hstack((distances[j], p_f_dists[j]))

        # for individual plots:
        # for j in range(len(densities)):
        #     genplots.plot_inv_cdf(
        #         axs[j], distances[j], yscale='log',
        #         label=labels[i], color=colors[i])
        group_distances.append(distances)

    indices = range(len(densities))  # [6,7,10] #density indices
    print densities
    for k, j in enumerate(indices):
        fig = plt.figure(figsize=(4, 3))
        p = densities[j]
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel(settings.get_prop_tex_name(settings.link_distance_tag))
        ax.set_ylabel(r"1-CDF(d)")
        print p, j
        ax.text(.5, 0.10, r"$\rho$ = " + str(p*100) + "\%",
                ha='center', va='center', transform=ax.transAxes)
        for i in range(len(fname_group_list)):
            genplots.plot_inv_cdf(
                ax, group_distances[i][j], label=labels[i],
                color=colors[i], yscale='log')
        plt.tight_layout()
        fig.savefig(cfg['outdata_dir'] +
                    "linkDistProbs_"+str(p)+".pdf", format="pdf", bbox_inches="tight")
        plt.close(fig)
    return None


def plot_pooled_corr_t_val_dists(cfg):
    """
    Plot the tvalue distributions for movie and rest

    Parameters
    ----------
    config : dict
        The following keys are required::

            "group_1_mat_fnames"
            "group_2_mat_fnames"
            "group_1_color"
            "group_2_color"
            "group_1_label"
            "group_2_label"
            "outdata_dir"
            "paired"

    Returns
    -------
    fig : the matplotlib figure object
    """
    config.require(cfg, ["group_1_mat_fnames", "group_2_mat_fnames",
                         "group_1_color", "group_2_color",
                         "group_1_label", "group_2_label",
                         "outdata_dir", "paired"]
                   )
    # get tvals
    flat_mats = dataio.get_blacklist_filtered_and_flattened_adj_mats(
        cfg["all_fnames"], cfg["blacklist_fname"]
    )
    if cfg["paired"]:
        t_vals = measures.paired_t_value(
            flat_mats, len(cfg['group_1_mat_fnames']))
    else:
        t_vals = measures.unpaired_t_value(
            flat_mats, len(cfg['group_1_mat_fnames']))

    minVal = np.min(t_vals)
    maxVal = np.max(t_vals)
    n_bins = 100
    bins, binCenters = aux.get_lin_bins(
        n_bins, minVal - (np.abs(minVal) * 0.1), maxVal +
        (np.abs(maxVal) * 0.1))
    bin_counts = aux.get_bin_counts(t_vals, bins)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # normalize
    bin_counts = bin_counts * 1. / (np.sum(bin_counts) * (bins[1] - bins[0]))
    labels = cfg['group_1_label'], cfg["group_2_label"]
    ax.plot(binCenters, bin_counts, label=r"" + labels[0] + "-" + labels[1])
    ax.set_xlabel(settings.get_prop_tex_name(settings.tval_tag))
    ax.set_ylabel(r"Probability density")
    ax.legend(loc=0)
    plt.savefig(cfg["outdata_dir"] + "tvalDist.pdf",
                format="pdf", bbox_inches='tight')
    return fig


def plot_global_w_props(cfg):
    """
    Plots the global network properties.

    Returns
    -------
    fig : a matplotlib figure
    """
    data = _get_global_props_data_to_plot(
        cfg, True, cfg['global_w_props'], True)
    colors = [cfg["group_1_color"], cfg["group_2_color"]]
    labels = [cfg["group_1_label"], cfg["group_2_label"]]
    return _global_props_plot(
        cfg["outdata_dir"],
        data,
        "globalWProps.pdf",
        colors,
        labels,
        paired=cfg['paired'],
        n=len(cfg["group_1_mat_fnames"]))


def plot_global_uw_props(cfg):
    """
    Plots the global network properties.

    Returns
    -------
    fig : a matplotlib figure
    """
    data = _get_global_props_data_to_plot(
        cfg, False, cfg['global_uw_props'], True)
    colors = [cfg["group_1_color"], cfg["group_2_color"]]
    labels = [cfg["group_1_label"], cfg["group_2_label"]]
    return _global_props_plot(
        cfg["outdata_dir"], data, "globalUWProps.pdf", colors,
        labels, paired=cfg['paired'], n=len(cfg["group_1_mat_fnames"]))


def _global_props_plot(
        outdata_dir,
        data,
        filenamePostFix,
        colors,
        labels,
        x_len=None,
        y_len=None, paired=False, n=-1):
    propMeans, propMeanConfIntervals, propPVals, densities = data

    x_len, y_len = _get_subplot_x_y(len(propMeans))
    fig = plt.figure(figsize=(x_len * 4, y_len * 4))
    color0 = colors[0]
    color1 = colors[1]
    for i, key in enumerate(propMeans):
        ax = host_subplot(y_len, x_len, i + 1)
        means0 = propMeans[key][0]
        means1 = propMeans[key][1]
        confInt0 = propMeanConfIntervals[key][0].T
        confInt1 = propMeanConfIntervals[key][1].T
        _plot_comparison_and_p_value(ax, densities, means0, means1, confInt0,
                                     confInt1, color0, color1,
                                     propPVals[key], key, paired=paired, n=n,
                                     labels=labels, xscale='log')
    plt.tight_layout(pad=2.2)
    fig.savefig(outdata_dir + filenamePostFix,
                format="pdf", bbox_inches='tight', pad_inches=0.5)
    return fig


def _plot_comparison_and_p_value(
        ax, densities, means0, means1, confInt0, confInt1, color0, color1,
        pVals=None, measure_key="", paired=False, n=-1, labels=None,
        xscale="linear", legendloc=0, significanceTresholdLine=None):

    ax.set_xlabel(settings.get_prop_tex_name(settings.densities_tag))
    ax.set_ylabel(settings.get_prop_tex_name(measure_key))
    # pax.set_ylabel("P-value")
    p1 = ax.plot(densities, means0, label=labels[0], color=color0, lw=1.5)
    p2 = ax.plot(densities, means1, label=labels[1], color=color1, lw=1.5)
    alpha = 0.25
    ax.fill_between(
        densities, confInt0[0], confInt0[1], color=color0, alpha=alpha)
    ax.fill_between(
        densities, confInt1[0], confInt1[1], color=color1, alpha=alpha)

    pax = ax.twinx()
    if significanceTresholdLine is not None:
        pax.axhline(significanceTresholdLine, xmin=0, xmax=1,
                    color="0.4", lw=1.0, ls="--", zorder=-1000)
    if paired:
        pax.axhline(y=2 ** -(n - 1), xmin=0, xmax=1,
                    color='0.6', lw=1.5, ls="-", zorder=-1000)
    p3 = pax.semilogy(densities, pVals, "-o", color="0.25",
                      markersize=2.5, label=r"$p$-value", zorder=-10)
    # to take the left side log scale off (a bug in matplotlib with loglog)
    pax.yaxis.tick_right()
    pax.set_ylim((10 ** -4, 10 ** 0))

    pax.set_xscale(xscale)
    ax.set_xscale(xscale)
    if xscale == "log":
        ax.set_xlim((np.min(densities), np.max(densities)))
        pax.set_xlim((np.min(densities), np.max(densities)))
    # make the legend
    lns = p1 + p2 + p3
    labs = [l.get_label() for l in lns]
    l = pax.legend(lns, labs, loc=legendloc, numpoints=1, handlelength=1)
    l.get_frame().set_alpha(0.5)
    return ax, pax


# def _temporalGlobalUWPropsPlot(fnGroups, bs_samples, bs_coverage, labels,
#                                props=settings.global_uw_props):
#     """
#     Plots the global uw properties with time. Each group in fnGroup denotes
#     a time point, in order of appearance.
#     Labels correspond to the times.
#     Parameters
#     ----------
#     fnGroups: groups of filenames representing the group at one time interval
#     bs_samples: number of bootstrap samples
#     bs_coverage: e.g. 95 (an integer) for 95 percent bootstrap coverage
#     labels: labels for the filenamegroups
#     props: properties to plot
#     outFileNamePrefix
#     """
#     propMeans, propMeanConfIntervals, _, densities = \
#         _get_global_props_data_to_plot(fnGroups, bs_samples, bs_coverage, False,
#                                  props, loadStatistics=False)
#     cmap = cm.get_cmap('jet')
# ["r", "b", "g", "purple", "yellow", "k"]
#     colors = cmap([np.arange(len(densities)) * 1. / len(densities)])[0]
#     x_len = 3
#     y_len = 2
#     for i, p in enumerate(densities):
#         fig = plt.figure(figsize=(10, 6))
#         fig.suptitle(settings.get_prop_tex_name(
#             settings.densities_tag) + str(p) + " globuwprops")
#         for j, key in enumerate(propMeans.keys()):
#             ax = fig.add_subplot(y_len, x_len, j)
#             _plot_mean_with_shaded_errors(ax, propMeans[key][:, i],
#                                       propMeanConfIntervals[key][:, i, :].T,
#                                       xticklabels=labels,
#                                       color=colors[i], xlabelRotation=45)
#             ax.set_ylabel(settings.get_prop_tex_name(key))
#             ax.grid()
#         plt.tight_layout(2.0)
#         fig.savefig(cfg["outdata_dir"] +
#                     "globUWProps_totnparts_" + str(len(fnGroups)) + "_density_"
#                     + str(p) + ".pdf")
#     fig = plt.figure(figsize=(15, 9))
#     suptitle = "globuwprops, densities:"
#     for p in densities:
#         suptitle = suptitle + " " + str(p)
#     fig.suptitle(suptitle)
#     for j, key in enumerate(propMeans.keys()):
#         ax = fig.add_subplot(y_len, x_len, j)
#         for i, p in enumerate(densities):
#             _plot_mean_with_shaded_errors(ax, propMeans[key][:, i],
#                                       propMeanConfIntervals[key][:, i, :].T,
#                                       xticklabels=labels, color=colors[i],
#                                       xlabelRotation=45)
#             ax.set_ylabel(settings.get_prop_tex_name(key))
#             ax.grid(True)
#         plt.tight_layout(2.0)
#     fig.savefig(cfg["outdata_dir"] +
#                 "globUWProps_totnparts_" + str(len(fnGroups)) + ".pdf")

def plotSameLinksShareVsLouvainSimilarityMeasures(cfg,
                                                  filenamesGroup1,
                                                  filenamesGroup2,
                                                  density=None):
    allFNames = filenamesGroup1 + filenamesGroup2
    data = dataio.mergeAndLoadLouvainProperties(allFNames, density)
    clusterings = data[settings.louvain_cluster_tag]
    simMatricesDict = gencomps.computeClusterSimilarityMeasures(clusterings)

    linkSimMatData = dataio.load_pickle(
        fnc.get_fname(cfg, settings.common_links_tag))
    index = 0
    densities = linkSimMatData[settings.densities_tag]
    if density is None:
        pass
    else:
        for j in range(len(densities)):
            if int(densities[index]) == int(density):
                index = j
                break
    linkSimMatrix = linkSimMatData[settings.common_links_tag][index]
    linkSimMatrix = linkSimMatrix / float(linkSimMatrix[0, 0])  # get to ratio

    triu_indices = np.triu_indices_from(linkSimMatrix, 1)
    for measure, simMatrix in simMatricesDict.iteritems():
        simMeasures = simMatrix[triu_indices]
        linkSims = linkSimMatrix[triu_indices]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(linkSims, simMeasures, "o")
        ax.set_ylabel(settings.get_prop_tex_name(measure))
        ax.set_xlabel(settings.get_prop_tex_name(settings.common_links_tag))
        (r, p) = pearsonr(linkSims, simMeasures)
        fig.suptitle(r"Pearson correlation: {:.3f}".format(r))
        fig.savefig(cfg["outdata_dir"] +
                    settings.common_links_tag +
                    "_vs_" + measure + ".pdf",
                    format="pdf", bbox_inches='tight')


def plotDensityVsAvgCommonFractionOfLinksOverAllPairs(cfg):
    data = dataio.get_fname(cfg, settings.common_links_tag)
    linkSimMatrices = data[settings.common_links_tag]
    densities = data[settings.densities_tag]
    linkSimAvgs = []
    linkSimStds = []
    triu_indices = np.triu_indices_from(linkSimMatrices[0], 1)
    for linkSimMatrix in linkSimMatrices:
        linkSimMatrix = linkSimMatrix / \
            float(linkSimMatrix[0, 0])  # get to ratio
        linkSimAvgs.append(np.average(linkSimMatrix[triu_indices]))
        linkSimStds.append(np.std(linkSimMatrix[triu_indices]))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    densities = np.array(densities)
    ax.errorbar(
        densities, linkSimAvgs, linkSimStds, label="avg. over all pairs")
    ax.set_xscale('log')
    l = ax.legend(loc=1)
    l.get_frame().set_alpha(0.75)
    if cfg['include_mst']:
        plt.savefig(cfg["outdata_dir"] + settings.common_links_tag
                    + "_vs_" + settings.densities_tag + "_with_mst.pdf",
                    format="pdf",
                    bbox_inches="tight")
    else:
        plt.savefig(cfg["outdata_dir"] + settings.common_links_tag
                    + "_vs_" + settings.densities_tag + ".pdf",
                    format="pdf", bbox_inches="tight")


def plotDensityVsAvgCommonFractionOfLinksOverAllPairsInDifferentGroups(
        cfg):
    data = dataio.load_pickle(fnc.get_fname(cfg, settings.common_links_tag))
    linkSimMatrices = data[settings.common_links_tag]
    densities = data[settings.densities_tag]
    linkSimAvgs = []
    linkSimStds = []
    n = len(linkSimMatrices[0]) / 2
    differentGroupIndices = (
        np.array((np.ones((n, n)) - np.eye(n)).nonzero()).T +
        np.array([n, n])).T
    print differentGroupIndices

    for linkSimMatrix in linkSimMatrices:
        linkSimMatrix = linkSimMatrix / \
            float(linkSimMatrix[0, 0])  # get to ratio
        linkSimAvgs.append(np.average(linkSimMatrix[differentGroupIndices]))
        linkSimStds.append(np.std(linkSimMatrix[differentGroupIndices]))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    densities = np.array(densities)
    ax.errorbar(densities, linkSimAvgs, linkSimStds,
                label="avg. over all pairs in different groups, neglecting" +
                " same subjects, stdev")
    ax.set_xscale('log')
    l = ax.legend(loc=1)
    l.get_frame().set_alpha(0.75)
    if cfg['include_mst']:
        plt.savefig(cfg["outdata_dir"] + settings.common_links_tag
                    + "_vs_" + settings.densities_tag +
                    "_between_groups.pdf",
                    format="pdf", bbox_inches="tight")
    else:
        plt.savefig(cfg["outdata_dir"] + settings.common_links_tag
                    + "_vs_" + settings.densities_tag +
                    "_between_groups_noMST.pdf",
                    format="pdf", bbox_inches="tight")


def plotDensityVsAvgPairedShare(cfg):
    """
    Plots the avg. fraction of the common links between the two settings
    (paired setup)
    """
    data = dataio.load_pickle(fnc.get_fname(cfg, settings.common_links_tag))
    linkSimMatrices = data[settings.common_links_tag]
    n = len(linkSimMatrices[0]) / 2
    densities = data[settings.densities_tag]
    linkSimAvgs = []
    linkSimStds = []
    vals = []
    for linkSimMatrix in linkSimMatrices:
        linkSimMatrix = linkSimMatrix / \
            float(linkSimMatrix[0, 0])  # get to ratio
        linkSimMatrix[range(n)]
        vals.append(linkSimMatrix[range(n), n + np.array(range(n))])
        linkSimAvgs.append(np.average(vals[-1]))
        linkSimStds.append(np.std(vals[-1]))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    densities = np.array(densities)
    ax.set_xlabel(settings.get_prop_tex_name(settings.densities_tag))
    ax.set_ylabel(r"Share of common links")
    vals = np.array(vals).T
    for i, val in enumerate(vals):
        ax.plot(densities, val, "r", alpha=0.4, label=str(i + 1))

    ax.errorbar(densities, linkSimAvgs, linkSimStds, color="b", label="avg")
    ax.set_xscale('log')
    l = ax.legend(loc=1)
    l.get_frame().set_alpha(0.75)
    if cfg['include_mst']:
        plt.savefig(cfg["outdata_dir"] + settings.common_links_tag
                    + "_paired_vs_" + settings.densities_tag + "_with_mst.pdf",
                    format="pdf", bbox_inches="tight")
    else:
        plt.savefig(cfg["outdata_dir"] + settings.common_links_tag
                    + "_paired_vs_" + settings.densities_tag + ".pdf",
                    format="pdf", bbox_inches="tight")


def plotLouvainSimilarityMeasures(cfg, additional_statistics=False,
                                  transparent=False):
    """
    Plots results for the (Louvain) community detection.

    Prints out also results for module size and Q differences into file
    outputDir/modularityStats.txt
    """
    n1 = len(cfg['group_1_mat_fnames'])
    n2 = len(cfg['group_2_mat_fnames'])
    data = dataio.mergeAndLoadLouvainProperties(
        cfg['group_1_mat_fnames'] + cfg['group_2_mat_fnames'], cfg['density'])
    clusterings = data[settings.louvain_cluster_tag]
    cluster_n = []
    for i, c in enumerate(clusterings):
        cluster_n.append(len(np.unique(c)))

    string = ""
    if additional_statistics:
        results_sizes = ptests.mean_difference_permtest(
            np.array(cluster_n), n1, n2, cfg['paired'],
            cfg['n_it_permutation']
        )
        results_mod = ptests.mean_difference_permtest(
            data[settings.modularity_tag], n1, n2, cfg['paired'],
            cfg['n_it_permutation'])
        string = "Cluster sizes:\n" + str(cluster_n) + "\n"
        string += str(np.ravel(results_sizes)) + "\n"
        string += "Cluster modularities: \n" + \
            str(np.ravel(data[settings.modularity_tag])) + "\n"
        string += str(np.ravel(results_mod)) + "\n"

    simMatricesDict = gencomps.computeClusterSimilarityMeasures(clusterings)

    for measure in simMatricesDict:
        if additional_statistics:
            string += "\n\n across group stats for measure"
            string += measure + \
                ": (crossMean, semidiagMean, semidiagMean - crossMean) \n"

            string += str(ptests.sim_matrix_semidiag_vs_inter_group_permtest(
                simMatricesDict[measure], 1e5, seed=10))
        fig = plt.figure(figsize=(8, 10))
        ax = fig.add_subplot(1, 1, 1)
        plot_sim_mat(ax, simMatricesDict[measure], measure)
        labels = [cfg['group_1_label'], cfg['group_2_label']]
        if ((cfg['paired'] is not False) and
                (cfg['n_it_permutation'] is not None) and
                (np.minimum(n1, n2) > 0) and
                (cfg['bootstrap_samples'] is not None) and
                (cfg['bootstrap-coverage'] is not None)):

            _addPermutationTestAndTextsToSimMatrix(
                fig, simMatricesDict[measure],
                measure, labels, n1, n2,
                cfg['paired'], cfg['n_it_permutation'], cfg[
                    'bootstrap_samples'],
                cfg['bootstrap_coverage'])
        ax.set_xticklabels(
            [str(i) for i in range(1, 1 + len(cfg['group_1_mat_fnames']))] +
            [str(i) for i in range(1, 1 + len(cfg['group_2_mat_fnames']))]
        )
        ax.set_yticklabels(
            [str(i) for i in range(1, 1 + len(cfg['group_1_mat_fnames']))] +
            [str(i) for i in range(1, 1 + len(cfg['group_2_mat_fnames']))]
        )
        fig.savefig(cfg["outdata_dir"] + "louvain_" +
                    measure + "_density_" + str(cfg['density']) + ".pdf", format="pdf", bbox_inches='tight')

    if additional_statistics:
        f = open(cfg["outdata_dir"] + "modularityStats.txt", "w")
        f.write(string)
        f.close()


def _addPermutationTestAndTextsToSimMatrix(
        fig, simMatrix, measure, labels, n1, n2,
        paired, permSamples, bs_samples, bs_coverage):
    commonX = 0.1
    commonY = 0.93
    controlX = 0.69
    controlY = 0.48
    treatmentX = 0.34
    treatmentY = 0.76
    fig.text(treatmentX, commonY, labels[0],
             rotation='horizontal', va='center', ha='center')
    fig.text(commonX, treatmentY, labels[
             0], rotation='vertical', va='center', ha='center')
    fig.text(controlX, commonY, labels[
             1], rotation='horizontal', va='center', ha='center')
    fig.text(commonX, controlY, labels[
             1], rotation='vertical', va='center', ha='center')
    print "permuation test statistics for the similarity measure " + measure
    stats = ptests.sim_matrix_within_group_mean_diff_permtests(
        simMatrix, n1, n2, paired, permSamples)
    print stats
    means, meanConfIntervals = \
        bootstrap.mean_groupwise_conf_intervals_from_sim_matrix(simMatrix, n1,
                                                                bs_samples,
                                                                bs_coverage)

    treatmentText = \
        (r"\begin!center?mean = {:.3f} \\ 95\%: ({:.3f}-{:.3f})\end!center?"
         ).format(means[0], meanConfIntervals[0][0], meanConfIntervals[0][1])
    controlText = \
        (r"\begin!center?mean = {:.3f} \\ 95\%: ({:.3f}-{:.3f})\end!center?"
         ).format(means[1], meanConfIntervals[1][0], meanConfIntervals[1][1])
    treatmentText = treatmentText.replace("!", "{")
    controlText = controlText.replace("!", "{")
    treatmentText = treatmentText.replace("?", "}")
    controlText = controlText.replace("?", "}")
    fig.text(treatmentX - 0.05, treatmentY - 0.05, treatmentText,
             bbox={"facecolor": "w", "alpha": 0.8}, va='center',
             ha='center', alpha=0.8)
    fig.text(controlX - 0.05, controlY - 0.05, controlText,
             bbox={"facecolor": "w", "alpha": 0.5}, va='center',
             ha='center')
    fig.text(0.50, 0.30, r"permutation test: p = {:.5f}".format(
        stats[0]), va='center', ha='center')


def plot_link_sim_mat(cfg):
    linkSimMatrix = dataio.load_pickle(
        fnc.get_fname(cfg, settings.common_links_tag))
    linkSimMatrix = linkSimMatrix / float(linkSimMatrix[0, 0])  # get to ratio
    fig = plt.figure(figsize=(8, 10))
    ax = fig.add_subplot(1, 1, 1)
    measure = settings.common_links_tag
    plot_sim_mat(ax, linkSimMatrix, measure)
    _addPermutationTestAndTextsToSimMatrix(fig, linkSimMatrix, measure)
    ax.set_xticklabels(2 * [str(i) for i in range(1, 14)])
    ax.set_yticklabels(2 * [str(i) for i in range(1, 14)])
    fig.savefig(cfg["outdata_dir"] + measure + ".pdf",
                format="pdf", bbox_inches='tight')


def _get_subplot_x_y(n):
    """ Computes somewhat appropriate x, y for fitting in n plots (x*y >= n)
    """
    startval = np.sqrt(n)
    y = np.floor(startval)
    x = np.ceil(n / y)
    assert x * y >= n
    return x, y


def plot_sim_mat(ax, m, measure, vmin=None, vmax=None):
    """
    Plot a similarity matrix in to the given axis "ax". "m" contains the image
    and "measure" is the name of the measure to be plotted (ie. a string)
    """
    sortedvals = np.sort(np.unique(m.flatten()))
    if measure in ['vi']:
        cmap = cm.hot
        if vmax is None:
            vmax = np.max(m)
        if vmin is None:
            try:
                vmin = sortedvals[1]
            except:
                vmin = sortedvals[0]
    elif measure in ['nmi', 'adjusted_rand', settings.common_links_tag]:
        cmap = cm.hot_r
        sortedvals = np.sort(np.unique(m.flatten()))
        if vmax is None:
            vmax_start = sortedvals[-1]
            for i in range(2, len(sortedvals)):
                vmax = sortedvals[-i]
                if (vmax_start - vmax) / vmax_start > 0.0001:
                    break
        if vmin is None:
            vmin = np.min(m)
    else:
        print "trying to plot unknown measure in function plot_sim_mat..." + \
            "defaulting the colormap hot_r"
        cmap = cm.hot_r
        if vmin is None:
            vmin = np.min(m)
        if vmax is None:
            vmax = np.max(m)

    im = ax.imshow(m, interpolation='nearest', cmap=cmap, vmax=vmax, vmin=vmin)
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal')
    cbar.set_label(settings.get_prop_tex_name(measure))

    ticks = np.arange(0, len(m))
    ax.yaxis.set_ticklabels(ticks + 1)
    ax.yaxis.set_ticks(ticks)
    ax.xaxis.set_ticks(ticks)
    ax.xaxis.set_ticklabels(ticks + 1)

    ax.xaxis.set_ticks_position('top')
    ax.tick_params(length=0, width=0, colors='k')
    ymed = np.average(ax.get_ylim())
    xmed = np.average(ax.get_xlim())
    ax.axhline(y=ymed, xmin=0, xmax=1, color='0.0', lw=1.0, ls="-")
    ax.axvline(x=xmed, ymin=0, ymax=1, color="0.0", lw=1.0, ls="-")


def _get_global_props_data_to_plot(cfg, weighted, props, load_stats=False):
    """
    Internal helper function, for fetching data for global
    network properties.

    Returns
    -------
    propMeans : a python dictionary with keys as the property
        propMeans[prop][mode][p] = mean
    propMeanConfIntervals : dict
        propData[prop][mode][p] = [low, high]
    propStats : a python dict containing the perm. ttest pvalues
        propStats[prop][p] = pvalue
    densities : list
        of the densities
    """
    if weighted:
        props_tag = 'global_w_props'
    else:
        props_tag = 'global_uw_props'

    fnGroups = [cfg['group_1_mat_fnames'], cfg['group_2_mat_fnames']]
    data = [dataio.merge_and_load_props_data(fnGroup, props_tag, props, cfg)
            for fnGroup in fnGroups]
    if load_stats:
        statsData = dataio.load_pickle(fnc.get_stats_fname(cfg, props_tag))

    propMeans = {}
    propMeanConfIntervals = {}
    propPVals = {}
    densities = data[0][settings.densities_tag]

    for prop in data[0].keys():
        if prop in props:
            modeMeans = []
            modeMeanErrs = []
            for i, datum in enumerate(data):
                # if prop != settings.densities_tag:
                propData = datum[prop]
                modeMeans.append(np.average(propData, 0))
                modeMeanErrs.append(
                    bootstrap.mean_conf_interval(
                        propData.T,
                        cfg['bootstrap_samples'],
                        cfg['bootstrap_coverage']
                    )
                )
            propMeans[prop] = np.array(modeMeans)
            propMeanConfIntervals[prop] = np.array(modeMeanErrs)
            if load_stats:
                propPVals[prop] = statsData[prop][settings.pval_tag]
    return propMeans, propMeanConfIntervals, propPVals, densities
