import dataio
import settings
import statistics
import numpy as np
import brainnets_params as params

def computeNodePVals():
    nodeProps = dataio.mergeAndLoadNodeProperties()
    nodePValStats = computePermutationTestStats(nodeProps, statistics.tValue, nIt=-1, paired='all')
    fName = dataio.getNodePValStatsFilename()
    dataio.savePickle(fName, nodePValStats)
    
def computeNodePropFDRStats(nIt = int(1e5), paired='all'):
    #load concatenated nodeprops
    nodeProps = dataio.mergeAndLoadNodeProperties() #1
    outDict = computeFDRStats(nodeProps, nIt, paired=paired)
    dataio.savePickle(dataio.getNodeFDRStatsFilename(), outDict)

def computeCorrFDRStatsPreDefinedThresholds():
    """
    The thresholds are fixed before. The data from this function is bound to
    be used for comparison with the mafdr function of matlab.
    
    Thresholds are symmetric.
    This function should be run on triton with 12 cores.
    (pFDR+estimation of null distribution)
    """
    ### load data
    corrData = dataio.loadFlattenedCorrMatrices()
    #fisher transformation: (in place to save memory!)
    np.arctanh(corrData, corrData)
    thresholds = np.linspace(2.,9.,15.)
    fdrs, _, _, ts = statistics.parallelPermFDREstimate(corrData, params.n1, thresholds, 
                                                        seed=123456, nCpus = params.parallel_cpus, paired='all')
    outputDict = {}
    outputDict[settings.fdrs_abbr] = fdrs
    outputDict[settings.thresholds_abbr] = ts
    outputDict[settings.tval_abbr] = statistics.tValue(corrData, params.n1)
    dataio.savePickle(settings.mrResultsDir+"corrFDRStatsPredefinedThresholds.pkl", outputDict)
    return
    
def computeCorrPVals():
    """
    Compute the permutation test based p-values for each link/correlation.
    be used for comparison with the mafdr function of matlab.
    
    Thresholds are symmetric.
    This function should be run on triton with 12 cores.
    (pFDR+estimation of null distribution)
    """
    ### load data
    corrData = dataio.mrLoadFlattenedCorrMatrices()
    #fisher transformation: (in place to save memory!)
    np.arctanh(corrData, corrData)
    pVals, tVals = statistics.parallelPermutationTest(statistics.tValue, corrData, params.n1, nCpus = params.parallel_cpus, paired='all')

    outputDict = {}
    outputDict[settings.pval_abbr] = pVals
    outputDict[settings.tval_abbr] = tVals
    dataio.savePickle(settings.mrResultsDir+"corrPValuesPredefinedThresholds.pkl", outputDict)
    return
    
#def computeCorrFDRStats():
#    ### load data
#    corrData = dataio.loadFlattenedCorrMatrices()
#    paired='all'
#    #fisher transformation: (in place to save memory!)
#    np.arctanh(corrData, corrData)
#    percents = np.logspace(-1, -6,6)*100.0 #percents, logspace base is 10
#    ts = statistics.tValueThresholds(corrData, params.n1, percents = percents, asym=False)
#    fdrs, _, _, ts = statistics.parallelPermFDREstimate(corrData, params.n1, ts, 
#                                                        seed=123456, nIt = nIt, 
#                                                        nCpus = params.parallel_cpus, 
#                                                        paired=paired, 
#                                                        asym=asym)
#    outputDict = {}
#    outputDict[settings.fdrs_abbr] = fdrs
#    outputDict[settings.thresholds_abbr] = ts
#    outputDict[settings.tval_abbr] = statistics.tValue(corrData, params.n1)
#    dataio.savePickle(dataio.getCorrFDRStatsFilename(paired=paired, asym=asym), outputDict)
    
    
def computeFDRStats(propsData, nCpus=params.parallel_cpus):
    """
    Compute FDR stats using the SAM-like method for estimating the FDR level
    for some thresholds.
    
    Works for looping over different properties. Used (only?) when computing
    node level stats.
    """    
    print "something might still be wrong with this computeFDRStats function - be cautious!"    
    outputDict = {}
    for key, dataArray in propsData.iteritems():
        #loop over different properties, but ignore the percentages prop
        if key == settings.percentages_abbr:
            outputDict[settings.percentages_abbr] = dataArray
        else:
            print "computing stats for prop " + key + "..."
            propDict = {}
#            percents = np.logspace(-1, -6,6)*100.0 #percents, logspace base is 10
#            ts = statistics.tValueThresholds(dataArray, settings.mrNsubjs, percents = percents, asym=asym)
            ts = np.linspace(1.,8.,15.)

            fdrs, _, _, ts = statistics.parallelPermFDREstimate(dataArray, params.n1,
                                                                ts, seed=123456, nIt = params.permSamples, 
                                                                nCpus=nCpus, paired=params.studyType, asym=False)
            propDict[settings.tval_abbr] = statistics.tValue(dataArray, params.n1)
            propDict[settings.fdrs_abbr] = fdrs
            propDict[settings.thresholds_abbr] = ts
            outputDict[key] = propDict
    return  outputDict

def computePermutationTestStats(propsData, statFunc):
    outputDict = {}
    statAbbr = ""
    if statFunc == statistics.tValue:
        statAbbr = settings.tval_abbr
    if statFunc == statistics.meanDifference:
        statAbbr = settings.md_abbr
        
    for key, dataArray in propsData.iteritems():
        #loop over different properties, but ignore the percentages prop
        if key == settings.percentages_abbr:
            outputDict[settings.percentages_abbr] = dataArray
        else:
            print "computing stats for prop " + key + "..."
            propDict = {}
                        
                        
            pVals, testStats = statistics.permutationTest(statFunc, propsData[key], params.n1, params., paired=studyPaired())
            propDict[settings.pval_abbr] = pVals
            propDict[statAbbr] = testStats
            outputDict[key] = propDict
    return  outputDict
    
def studyPaired():
    if params.studyType == 'paired':
        return True
    return False


def computeGlobalWPropStats(statFunc=statistics.meanDifference):
    globalWProps = dataio.mrMergeAndLoadGlobalWProperties()
    outDict = computePermutationTestStats(globalWProps, statFunc)
    dataio.savePickle(settings.mrGetGlobalWStatsFilename(), outDict)
    
def computeGlobalUWPropStats(statFunc = statistics.meanDifference):
    globalUWProps = dataio.mrMergeAndLoadGlobalUWProperties()
    outDict = computePermutationTestStats(globalUWProps, statFunc)
        dataio.savePickle(settings.mrGetGlobalUWStatsFilename(), outDict)
    
    

