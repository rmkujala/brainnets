#
# General settings common for the usage of the module brainnets
# If something is added here, make sure the necessary requirements
# are done also in the program code.
# This module is practically imported by every other module of the
# library
#
# Authors: Rainer.Kujala@gmail.com
#

be_verbose = True

config_tag = "config"

# Statistics stuff:
density_tag = "density"
densities_tag = "density_range"
thresholds_tag = "thresh"
fdrs_tag = "fdrs"
pfdrs_tag = "pfdrs"
pi0_tag = "pi0"  # estimated proportion of the null distribution
qval_tag = "q-value"
tval_tag = "tvalues"
meandifference_tag = "mean_difference"
pval_tag = "pvalues"
mean_tag = "mean"
meanerror_tag = "meanerr"

# Link stuff
correlation_tag = "correlation"
correlation_distribution_tag = "corr_distr"
link_distance_tag = "link_distance"
link_distance_common_tag = "common_link_distances"
weight_tag = "weight"
avg_weight_tag = "avg_weight"
common_links_tag = "common_links"


# Global unweighted properties:
avg_clustering_tag = "average_clustering"
global_clustering_tag = "global_clustering"
avg_path_length_tag = "average_path_length"
efficiency_tag = "efficiency"
max_kshell_tag = "max_kshell"
assortativity_tag = "assortativity"
max_degree_tag = "max_degree"
global_uw_props = [avg_clustering_tag, global_clustering_tag,
                   avg_path_length_tag, max_kshell_tag,
                   assortativity_tag, max_degree_tag,
                   efficiency_tag]

# Global weighted properties
weighted_average_path_length_tag = "weighted_average_path_length"
max_strength_tag = "max_strength"
weighted_clustering_tag = "weighted_clustering"
weighted_assortativity_tag = "weighted_assortativity"
global_w_props = [weighted_average_path_length_tag, max_strength_tag,
                  weighted_clustering_tag, weighted_assortativity_tag,
                  avg_weight_tag]

# Node properties
degree_tag = "degree"
strength_tag = "strength"
betweenness_centrality_tag = "betweenness_centrality"
weighted_betweenness_centrality_tag = "weighted_betweenness_centrality"
k_shell_tag = "k_shell"
node_clustering_tag = "node_clustering"
node_props = [degree_tag, strength_tag, betweenness_centrality_tag,
              weighted_betweenness_centrality_tag, k_shell_tag,
              node_clustering_tag]
nodePropsWOMST = [degree_tag, strength_tag, k_shell_tag, node_clustering_tag]
# so far not in use:
# linkProps = ["edge_betweenness_centrality",
# "weighted_edge_betweenness_centrality"]

# Louvain modularity properties
modularity_tag = "modularity_q"
louvain_cluster_tag = "louvain_clusters"
louvain_consensus_tag = "louvain_consensus_clusters"
louvain_consensus_si_tag = "louvain_consensus_si"
louvainProps = [modularity_tag, louvain_cluster_tag]
undef_clu_label = -1

# Clustering comparisons:
vi_tag = "vi"
nmi_tag = "nmi"
adjusted_rand_tag = "adjusted_rand"

cluster_similarity_measures = ["vi", "nmi", "adjusted_rand"]
scaled_inclusivity_tag = "SI"


import brainnets
package_dir = brainnets.__path__[0] + "/"

_prop_tex_names = {
    avg_clustering_tag: r"Average clustering $C_a$",
    global_clustering_tag: r"Global clustering $C_g$",
    avg_path_length_tag: r"Avg. path length $\langle l \rangle$",
    assortativity_tag: r"Assortativity $r$",
    efficiency_tag: r"Efficiency $E$",
    max_degree_tag: r"Max. degree $k_{max}$",
    meandifference_tag: r"Mean difference",
    weighted_average_path_length_tag: r"Weighted avg. path length $l^w$",
    max_strength_tag: r"Max. strength $s_{max}$",
    max_kshell_tag: r"Max. k-shell $k^s_{max}$",
    weighted_clustering_tag: r"Weighted clustering $C_w$",
    weighted_assortativity_tag: r"Weighted assortativity $r^w$",
    degree_tag: r"Degree $k$",
    strength_tag: r"Strength $s$",
    betweenness_centrality_tag: r"Node betweenness $b_v$",
    weighted_betweenness_centrality_tag: r"Weighted node betweenness $b_v^w$",
    k_shell_tag: r"$k$-shell index $k^s$",
    weight_tag: r"Weight $w$",
    # edge_betweenness_centrality_tag: r"Edge betweenness $b_e$",
    # weighted_edge_betweenness_centrality":
    #           r"Weighted edge betweenness $b_e^w$",
    densities_tag: r"Network density $\rho$",
    modularity_tag: r"Modularity $Q$",
    vi_tag: r"Variation of information $VI$",
    nmi_tag: r"Normalized mutual information $NMI$",
    adjusted_rand_tag: r"Adjusted rand index $R$",
    correlation_tag: r"Correlation $c$",
    tval_tag: r"t-value $t$",
    common_links_tag: r"Fraction of common links",
    link_distance_tag: r"Distance d [mm]",
    avg_weight_tag: r"Average weight $\langle w \rangle$",
    node_clustering_tag: r"Local clustering coefficient $c_i$"
}

# directory where all created data is assumed to be located and outputted
# you can change this in your script if you wish


def get_prop_tex_name(key):
    """
    Obtain a tex_label for a key, if such is available.
    Otherwise returns the key itself.

    Parameters
    ----------
    key : str

    Returns
    -------
    tex_name: str
    """
    try:
        return _prop_tex_names[key]
    except:
        # just return the key if _prop_tex_names[key] was not found.
        return key
