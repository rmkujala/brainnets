
# a problem in triton sometimes?
try:
    from scipy.stats import ks_2samp
except:
    print "ks_2samp could not be imported"
import numpy as np

from brainnets import dataio, settings
from brainnets import fname_conventions as fnc
from verkko.permtests import measures
from verkko.permtests import ptests


def comp_node_pvals(cfg, matlabQValComps=False):
    """
    Computes the permutation-based p-values for the precomputed nodeproperties
    The comparison is made between
    ``cfg["group_1_mat_fnames"]`` and ``cfg["group_1_mat_fnames"]``


    Parameters
    ----------
    cfg : dict
        brainnets config dict
    """
    if not cfg['paired']:
        assert cfg['n_it_permutation'] != "all", \
            "all is not allowed for non-paired tests"

    fnames_group1 = cfg['group_1_mat_fnames']
    fnames_group2 = cfg['group_2_mat_fnames']
    fnames = fnames_group1 + fnames_group2

    node_props_data = dataio.merge_and_load_props_data(
        fnames, "node_props", cfg['node_props'], cfg)
    ok_nodes = dataio.get_ok_nodes(cfg['blacklist_fname'])
    for key in node_props_data:
        if key != settings.densities_tag:
            node_props_data[key] = node_props_data[key][:, ok_nodes]
    n1 = len(fnames_group1)
    n2 = len(fnames_group2)

    node_pval_stats = _comp_perm_test_stats_for_some_props(
        node_props_data, n1, n2, cfg['paired'], cfg['n_it_permutation'],
        measures.mean_difference, matlabQValComps=False)

    # print node_pval_stats.keys()
    # for nodeProp in settings.node_props:
    #     if nodeProp == settings.densities_tag:
    #         continue
    #     propPValStats = node_pval_stats[nodeProp]
    #     data = node_pval_stats[nodeProp]
    #     out_dict[settings.pval_tag] = \
    #         dataio.expand_1D_node_vals_to_non_blacklisted_array(
    #             propPValStats[settings.pval_tag], ok_nodes)
    #     out_dict[settings.meandifference_tag] = \
    #         dataio.expand_1D_node_vals_to_non_blacklisted_array(
    #             propPValStats[settings.meandifference_tag], ok_nodes)
    # if matlabQValComps:
    # out_dict[settings.pfdrs_tag] = \
    # dataio.expand_1D_node_vals_to_non_blacklisted_array(
    # propPValStats[settings.pfdrs_tag], ok_nodes)
    # out_dict[settings.qval_tag] = \
    # dataio.expand_1D_node_vals_to_non_blacklisted_array(
    # propPValStats[settings.qval_tag], ok_nodes)
    #     except:
    #         raise
    out_fname = fnc.get_stats_fname(cfg, "node_props")
    dataio.save_pickle(out_fname, node_pval_stats)


def comp_glob_prop_stats(cfg, weighted):
    """
    Computes the permutation-based p-values for the previously computed
    globalWPropStats

    Parameters
    ----------
    cfg : dict
        the brainnets config dictionary
    """
    if not cfg['paired']:
        assert cfg['n_it_permutation'] != "all", \
            "all is not allowed for non-paired tests"

    fnames_group1 = cfg['group_1_mat_fnames']
    fnames_group2 = cfg['group_2_mat_fnames']
    fnames = fnames_group1 + fnames_group2
    n1 = len(fnames_group1)
    n2 = len(fnames_group2)

    if weighted:
        props_tag = 'global_w_props'
    else:
        props_tag = 'global_uw_props'

    props = cfg[props_tag]
    global_props_data = dataio.merge_and_load_props_data(
        fnames, props_tag, props, cfg
    )
    out_dict = _comp_perm_test_stats_for_some_props(
        global_props_data, n1, n2, cfg['paired'], cfg['n_it_permutation'],
        measures.mean_difference, matlabQValComps=False)
    out_fname = fnc.get_stats_fname(cfg, props_tag)
    dataio.save_pickle(out_fname, out_dict)


def comp_simple_corr_dist_stats(cfg):
    """
    Compute correlation distribution statistics for the difference of the two
    modes.

    Parameters
    ----------
    cfg : cfg
        the brainnets config dict

    Returns
    -------
    ks_p : float
        the p-value for the result of two sided ks-test
    mean_difference : float
        the difference in the mean of the two distributions
    md_p :  float
        The two-sided p-value (based on permutation test) for the mean
        difference
    """
    if not cfg['paired']:
        assert cfg['n_it_permutation'] != "all", \
            "all is not allowed for non-paired tests"

    fnames_group1 = cfg['group_1_mat_fnames']
    fnames_group2 = cfg['group_2_mat_fnames']
    fnames = fnames_group1 + fnames_group2
    n1 = len(fnames_group1)
    n2 = len(fnames_group2)

    # load data
    corr_data = dataio.get_blacklist_filtered_and_flattened_adj_mats(
        fnames,
        cfg['blacklist_fname']
    )
    # load correlation stats:
    corr_data = dataio.get_blacklist_filtered_and_flattened_adj_mats(
        fnames_group1 + fnames_group2, cfg['blacklist_fname'])
    avgs = np.average(corr_data, axis=1)
    md_p, stats = ptests.mean_difference_permtest(
        avgs, n1, n2, cfg['paired'], cfg['n_it_permutation'], seed=123456)
    treatment_corrs = np.ravel(corr_data[:n1])
    control_corrs = np.ravel(corr_data[n1:])
    D, ks_p = ks_2samp(treatment_corrs, control_corrs)
    result_dict = {'D': D, 'ks_p': ks_p,
                   'mean_difference': stats, 'md_p': md_p}
    dataio.save_pickle(
        fnc.get_stats_fname(cfg, settings.correlation_distribution_tag),
        result_dict
    )


def comp_corr_p_vals(cfg):
    # fname_group1, fname_group2, blFName, paired, n_it, =1):
    """
    Compute the permutation test based p-values for each link/correlation.
    be used for comparison with the mafdr function of matlab.

    This function should be run on triton (with 12 cores).
    (pFDR+estimation of null distribution)
    """
    if not cfg['paired']:
        assert cfg['n_it_permutation'] != "all", \
            "all is not allowed for non-paired tests"

    fnames_group1 = cfg['group_1_mat_fnames']
    fnames_group2 = cfg['group_2_mat_fnames']
    fnames = fnames_group1 + fnames_group2
    n1 = len(fnames_group1)
    n2 = len(fnames_group2)

    # load data
    corr_data = dataio.get_blacklist_filtered_and_flattened_adj_mats(
        fnames,
        cfg['blacklist_fname']
    )
    # fisher transformation: (in place to save memory!)
    np.arctanh(corr_data, corr_data)
    print "starting permutation testing"
    pVals, tVals = ptests.mean_difference_permtest(
        corr_data, n1, n2, cfg['paired'], cfg['n_it_permutation'], seed=123456,
        n_cpus=cfg['n_cpus'])
    out_dict = {}
    out_dict[settings.pval_tag] = pVals
    out_dict[settings.tval_tag] = tVals
    dataio.save_pickle(fnc.get_stats_fname(cfg, settings.correlation_tag),out_dict)


# def computeCorrQVals():
#     """ Compute in hammer """
#     pVals = dataio.loadPickle(
#         dataio.getCorrPValStatsFileName())[settings.pval_tag]
#     fdrs, qvals, pi0 = statistics.pValuesToPFDRWithMatlab(pVals)
#     out_dict = {}
#     out_dict[settings.pfdrs_tag] = fdrs
#     out_dict[settings.qval_tag] = qvals
#     out_dict[settings.pi0_tag] = pi0
#     dataio.save_pickle(dataio.getCorrQValStatsFileName(), out_dict)


# permutation based FDR and pFDR computations:
#

# def computeNodePropFDRStats(fname_group1, fname_group2, blacklistFileName,
    # paired, n_it, nVoxels,  percentage=None, n_cpus=1, outFileName=None):
# load concatenated nodeprops
#     nodePropsData = dataio.mergeAndLoadNodeProperties(
#         fname_group1 + fname_group2, percentage)
#     n1 = len(fname_group1)
#     n2 = len(fname_group2)
#     ok_nodes = dataio.get_ok_nodes(nVoxels, blacklistFileName)
#     for nodeProp, data in nodePropsData.iteritems():
#         try:
#             nodePropsData[nodeProp] = data[:, ok_nodes]
#         except:
#             if nodeProp == settings.densities_tag:
#                 continue
#             raise
#     out_dict = computeFDRStats(nodePropsData, n1, n2, paired, n_it, n_cpus=1)
#     for prop, fdrdict in out_dict.iteritems():
#         if prop != settings.densities_tag:
#             out_dict[prop][settings.tval_tag] = \
#                   dataio.expand_1D_node_vals_to_non_blacklisted_array(
#                 fdrdict[settings.tval_tag], ok_nodes, default_value \
                        # =float("nan"))
#     if outFileName == None:
#         dataio.save_pickle(dataio.getNodeFDRStatsFileName(percentage),
                             # out_dict)
#     else:
#         dataio.save_pickle(outFileName, out_dict)


# def computeCorrFDRStats(fname_group1, fname_group2, blFName, paired, n_it,
                          # n_cpus=1):
#     """
#     The thresholds are fixed before. The data from this function is bound to
#     be used for comparison with the mafdr function of matlab.
#     """
# load data
# outFname = dataio.getCorrFDRStatsFileName()  # so the function works..
#     corr_data = dataio.loadFlattenedAdjacencyMatrices(
#         fname_group1 + fname_group2, blFName)
# fisher transformation: (in place to save memory!)
#     np.arctanh(corr_data, corr_data)
#     propDict = {settings.correlation_tag: corr_data}
#     n1 = len(fname_group1)
#     n2 = len(fname_group2)
#     out_dict = computeFDRStats(propDict, n1, n2, paired, n_it, n_cpus=n_cpus)
#     dataio.save_pickle(outFname, out_dict)


# def computeFDRStats(propsData, n1, n2, paired, n_it, n_cpus=1, asym=False,
#                     symThresholds=None):
#     """
#     Compute FDR stats using the SAM-like method for estimating the FDR level
#     for some thresholds.

#     Works for looping over different properties.
#     """
#     out_dict = {}
#     for key, dataArray in propsData.iteritems():
# loop over different properties, but ignore the percentages prop
#         if key == settings.densities_tag:
#             out_dict[settings.densities_tag] = dataArray
#         else:
#             print "computing stats for prop " + key + "..."
#             propDict = {}
# percents = np.logspace(-1, -6,6)*100.0 #percents, logspace base is 10
# ts = statistics.t_valueThresholds(dataArray, settings.mrNsubjs, percents = \
#        percents, asym=asym)
#             if symThresholds == None:
#                 symThresholds = np.linspace(2., 6., 11.)
# if asym:
# ts = np.array(list(ts)+list(-ts))
#             output = statistics.permTValueFDREstimate(
# dataArray, n1, n2, paired, n_it, symThresholds=symThresholds,
# asymThresholds=[], seed=1234, n_cpus=n_cpus)

#             [FDRsSym, pFDRsSym, FDRsAsym, pFDRsAsym] = output

#             if paired:
#                 propDict[settings.tval_tag] = statistics.pairedTValue(
#                     dataArray, n1)
#             else:
#                 propDict[settings.tval_tag] = \
#                     statistics.t_value(dataArray, n1)

#             propDict[settings.fdrs_tag] = FDRsSym
#             propDict[settings.pfdrs_tag] = pFDRsSym
#             propDict[settings.thresholds_tag] = symThresholds
#             out_dict[key] = propDict
#     return out_dict


def _comp_perm_test_stats_for_some_props(props_data, n1, n2, paired, n_it,
                                         stat_func, seed=10, n_cpus=1,
                                         matlabQValComps=False):
    out_dict = {}
    statTag = ""
    permutationTestFunc = None
    if stat_func in [measures.unpaired_t_value,
                     measures.paired_t_value]:
        statTag = settings.tval_tag
        permutationTestFunc = ptests.t_value_permtest
    if stat_func == measures.mean_difference:
        statTag = settings.meandifference_tag
        permutationTestFunc = ptests.mean_difference_permtest

    for key, dataArray in props_data.iteritems():
        # loop over different properties, but ignore the percentages prop
        if key == settings.densities_tag:
            out_dict[settings.densities_tag] = dataArray
        else:
            print "computing stats for prop " + key + "..."
            propDict = {}
            pVals, testStats = permutationTestFunc(
                props_data[key], n1, n2, paired, n_it, seed, n_cpus=n_cpus)
            propDict[settings.pval_tag] = pVals
            propDict[statTag] = testStats
            # if matlabQValComps:
            #     import statistics
            #     pfdrs, qvals, pi0 = statistics.pValuesToPFDRWithMatlab(pVals)
            #     propDict[settings.pfdrs_tag] = pfdrs
            #     propDict[settings.qval_tag] = qvals
            #     propDict[settings.pi0_tag] = pi0
            out_dict[key] = propDict
    return out_dict
