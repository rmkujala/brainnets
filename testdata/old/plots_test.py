from brainnets import settings, plots
import test_params as tp
import os
import numpy as np

clearAll = False#True

path = os.path.dirname(tp.__file__)+"/"
print path
settings.setBrainnetsDataDir(path)
settings.be_verbose = False

f = open(path+"plotTestLog.txt", "w")

#node properties
try:
    fnGroups = [tp.movieFNames, tp.restFNames]
    f.write("stateted plotting\n")
#    plots.plotLouvainSimilarityMeasures(tp.movieFNames, tp.restFNames, tp.density, tp.paired, tp.nIt, tp.labels, tp.bsSamples, tp.bsCoverage, transparent=False, additional_statistics=True)
#    plots.plotDensityVsAvgCommonFractionOfLinksOverAllPairs(tp.includeMST)
#    plots.plotDensityVsAvgCommonFractionOfLinksOverAllPairsInDifferentGroups(tp.includeMST)
#    plots.plotSameLinksShareVsLouvainSimilarityMeasures(tp.restFNames, tp.movieFNames)
    plots.globalWPropsPlot(tp.movieFNames, tp.restFNames,  tp. bsSamples, tp.bsCoverage, tp.labels, tp.colors)
    plots.globalUWPropsPlot(tp.movieFNames, tp.restFNames,  tp. bsSamples, tp.bsCoverage, tp.labels, tp.colors)
    plots.plotCorrTValDist(tp.movieFNames, tp.restFNames, tp.blFName, tp.labels, tp.colors, tp.paired)
    plots.plotCorrelationDists(fnGroups, tp.blFName, tp.labels, tp.colors)
    plots.plotIndividualCorrelationDists(fnGroups, tp.blFName, tp.colors)
    plots.plotLinkDistanceProbabilities(fnGroups, tp.labels, tp.colors)
    plots.plotCorrelationDistsSubjectwise(fnGroups, tp.blFName, tp.labels, tp.colors)
    f.write("finished plotting\n\n")
except:
    f.write("\n ...an error occurred...")
    f.close()
    raise


#rm all filenames with ending .pkl
if clearAll:
    import subprocess
    import glob
    files = glob.glob(path+"*.pkl")
    files+= glob.glob(path+"*.pdf")
    command = ['rm']+files
    #print command
    subprocess.call(command)
    
