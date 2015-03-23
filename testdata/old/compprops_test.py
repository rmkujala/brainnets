from brainnets import compprops
from brainnets import settings
from brainnets import compstats
import test_params as tp
import os
import numpy as np

clearAll = True

path = os.path.dirname(tp.__file__) + "/"
print path
settings.setBrainnetsDataDir(path)
settings.be_verbose = False

f = open(path + "compTestLog.txt", "w")

includeMST = tp.includeMST

n_cpus = 1  # 2 for parallellism
density = tp.density
nIt = tp.nIt

# node properties
try:
    f.write("computing nodeProperties\n")
    compprops.computeNodeProperties(
        tp.allFNames, tp.allBlFNames, tp.density, includeMST, n_cpus=n_cpus)
    compstats.comp_node_pvals(
        tp.movieFNames, tp.restFNames, tp.blFName, tp.nVoxels, True, nIt, density)
    f.write("finished computing nodeProperties\n\n")

    f.write("computing uw properties\n")
    compprops.computeGlobalUWProperties(
        tp.allFNames, tp.allBlFNames, tp.pRange, includeMST, n_cpus=n_cpus)
    compstats.computeGlobalUWPropStats(
        tp.movieFNames, tp.restFNames, True, 'all')
    f.write("computed uw properties\n\n")

    f.write("computing globalWeightedProperties\n")
    compprops.computeGlobalWProperties(
        tp.allFNames, tp.allBlFNames, tp.pRange, includeMST,  n_cpus=n_cpus)
    compstats.computeGlobalWPropStats(
        tp.movieFNames, tp.restFNames, False, nIt)
    f.write("finished globalWeightedProperties\n\n")

    f.write("computeLouvainCommunityProperties\n")
    compprops.computeLouvainCommunityProperties(
        tp.allFNames, tp.allBlFNames, tp.density, includeMST, nIt, n_cpus=n_cpus)
    f.write("finished computeLouvainCommunityProperties\n\n")

    f.write("computeLinkSimilarityMatrices\n")
    compprops.computeLinkSimilarityMatrices(
        tp.allFNames, tp.blFName, tp.pRange, includeMST)
    f.write("finished computeLinkSimilarityMatrices\n\n")

    f.write("computeLinkDistances \n")
    compprops.computeLinkDistances(
        tp.allFNames, tp.blFName, tp.pRange, tp.voxelInfoFileName, includeMST)
    f.write("finished computeLinkDistances \n\n")

    f.write("corrdist stats \n")
    compstats.comp_simple_corr_dist_stats(
        tp.movieFNames, tp.restFNames, tp.blFName, tp.paired, tp.nIt)
    f.write("finished corrdist stats \n\n")

    f.write("corr pvals \n")
    compstats.computeCorrPVals(
        tp.movieFNames, tp.restFNames, tp.blFName, tp.paired, tp.nIt, n_cpus)
    compstats.computeCorrQVals()
    f.write("finished corr pvals \n\n")

    f.write("nodePropFDRStats\n")
    compstats.computeNodePropFDRStats(
        tp.movieFNames, tp.restFNames, tp.blFName, tp.paired, tp.nIt, tp.nVoxels, tp.density, n_cpus=n_cpus)
    f.write("finished nodePropFDRStats\n\n")

    f.write("corrFDRStats\n")
    compstats.computeCorrFDRStats(
        tp.movieFNames, tp.restFNames, tp.blFName, tp.paired, tp.nIt, n_cpus=n_cpus)
    f.write("finished corrFDRStats\n")

    f.close()


except:
    f.write("\n ...an error occurred...")
    f.close()
    raise


# rm all filenames with ending .pkl
if clearAll:
    import subprocess
    import glob
    files = glob.glob(path + "*.pkl")
    command = ['rm'] + files
    # print command
    subprocess.call(command)
