#
# General settings common for the usage of the module brainnets
# If something is added here, make sure the necessary requirements 
# are done also in the program code.
# This module is practically imported by every other module of the
# library
#
# Authors: Rainer.Kujala@gmail.com
#

percentages_tag = "perc"
thresholds_tag = "thresh"
fdrs_tag = "fdrs"
tval_tag = "tvalues"
md_tag = "mean_difference"
pval_tag = "pvalues"
mean_tag = "mean"
meanerror_tag = "meanerr"
modularity_tag = "modularity_q"
louvain_cluster_tag = "louvain_clusters"
correlation_tag = "correlation"

#other tags:
globUWProps = ["average_clustering", "global_clustering", "average_path_length", "max_kshell", "assortativity", "max_degree"]
globWProps = ["weighted_average_path_length", "max_strength", "weighted_clustering", "weighted_assortativity"]
nodeProps = ["degree", "strength", "betweenness_centrality", "weighted_betweenness_centrality", "k_shell"]                
#so far not in use:
#linkProps = ["edge_betweenness_centrality", "weighted_edge_betweenness_centrality"]
louvainProps = [modularity_tag, louvain_cluster_tag]

cluster_similarity_measures = ["vi", "nmi", "adjusted_rand"]

_propTexNames = {
                "average_clustering": r"Average clustering $C_a$", 
                "global_clustering": r"Global clustering $C_g$", 
                "average_path_length": r"Avg. path length $\langle l \rangle$", 
                "assortativity": r"Assortativity $r$", 
                "max_degree": r"Max. degree $k_{max}$",
                "mean_difference": r"Mean difference",
                "weighted_average_path_length": r"Weighted avg. path length $l^w$", 
                "max_strength": r"Max. strength $s_{max}$", 
                "max_kshell": r"Max. k-shell $k^s_{max}$",
                "weighted_clustering": r"Weighted clustering $C_w$",
                "weighted_assortativity": r"Weighted assortativity $r^w$",
                "degree": r"Degree $k$", 
                "strength": r"Strength $s$", 
                "betweenness_centrality": r"Node betweenness $b_v$",
                "weighted_betweenness_centrality": r"Weighted node betweenness $b_v^w$",
                "k_shell": r"K-shell index $k^s$",
                "weight": r"Weight $w$",
                "edge_betweenness_centrality": r"Edge betweenness $b_e$",
                "weighted_edge_betweenness_centrality": r"Weighted edge betweenness $b_e^w$",
                "perc":r"Network density $\rho$ (in \%)",
                "q": r"Modularity $Q$",
                "vi": r"Variation of information $VI$",
                "nmi": r"Normalized mutual information $NMI$", 
                "adjusted_rand": r"Adjusted rand index $R$",
                "corr": r"Correlation $c$",
                "tvalues": r"t-value $t$"                
                }                

def getPropTexName(key):
    try:
        return _propTexNames[key]
    except:
        #return the key as a 
        return key
