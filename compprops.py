import multiprocessing
import dataio
import settings
import gencomps
import sys
import gc
import netgen
import brainnets_params as params

#The structure of the module could be improved, although it works.
#Everything is made so that the parallellism part was needed to code only once.
#Functions starting with an underscore are (at least originally) meant to be 
#used only internally.

#FUNCTIONS FOR THE END-USER:
def computeAndSaveNodeProperties(nCpus=params.parallel_cpus):
    print "starting computing node properties..."
    _loopOverModesAndSubjectsParallel(_nodePropWorkFunc, nCpus)

def computeAndSaveGlobalUWProperties(nCpus=params.parallel_cpus):
    _loopOverModesAndSubjectsParallel(_globalUWPropWorkFunc, nCpus)

def computeAndSaveGlobalWProperties(nCpus=params.parallel_cpus):
    _loopOverModesAndSubjectsParallel(_globalWPropWorkFunc, nCpus)
    
def computeAndSaveLouvainCommunityProperties(nCpus=params.parallel_cpus):
    _loopOverModesAndSubjectsParallel(_louvainClusteringPSpecialWorkFunc, nCpus)
    
def computeAndSaveGlobalWPropertiesTreatment(nCpus=params.parallel_cpus):
    """
    Compute global weighted properties for the treatment group.
    """
    _loopOverModesAndSubjectsParallel(_globalWPropWorkFunc, nCpus, [params.modes[0]])
    
def computeAndSaveGlobalWPropertiesControl(nCpus=params.parallel_cpus):
    """
    Compute global weighted properties for the control group.
    """
    _loopOverModesAndSubjectsParallel(_globalWPropWorkFunc, nCpus, [params.modes[1]])
    
#THE FUNCTION WHERE ACTUAL PARALLLELLISM TAKES PLACE:
def _loopOverModesAndSubjectsParallel(workFunc, nCpus, modes=None):
    """
    Loops over the two different modes in parallel using the given workFunc
    Args:
        workFunc: a *WorkFunc
            
    """
    inputs = []
    if modes == None: #enable looping over one mode only
        modes = params.modes
    for mode in modes:
        for i in range(1,1+params.modeToN[mode]):
            inputs.append( (mode, i) )
    pool = multiprocessing.Pool(processes=nCpus)
    # hack for keyboard interrupt, nothing should last more than a year 
    # (the digit is seconds)
    pool.map_async(workFunc, inputs, chunksize=1).get(31536000)

#"INTERNAL FUNCTIONS":
def _doStart(inData):
    """ Prints which work has started and returns the mode and i"""
    mode = inData[0]
    i = inData[1]
    print "started", mode,  i
    sys.stdout.flush()
    adjMat = dataio.loadFilteredAdjacencyMatrix(mode, i)
    return mode, i, adjMat
        
def _doEnd(mode, i, outFileName, results):
    """ Prints which works has ended and saves the results"""
    dataio.savePickle(outFileName, results)
    print "finished", mode,  i
    sys.stdout.flush()    

def _globalWPropWorkFunc(inData):
    mode, i, adjMat = _doStart(inData)
    results = gencomps.getGlobalPropsForPRange(adjMat, params.pRange, settings.globWProps, True)
    outFileName = dataio.getGlobalWPropsIndividualFileName(mode, i)
    _doEnd(mode, i, outFileName, results)

def _globalUWPropWorkFunc(inData):
    mode, i, adjMat = _doStart(inData)
    results = gencomps.getGlobalPropsForPRange(adjMat, params.pRange, settings.globUWProps, False)
    outFileName = dataio.getGlobalUWPropsIndividualFileName(mode, i)
    _doEnd(mode, i, outFileName, results)

def _nodePropWorkFunc(inData):
    mode, i, adjMat = _doStart(inData) 
    results = gencomps.getNodePropsForP(adjMat, params.pRange, props=settings.nodeProps)
    outFileName = dataio.getNodePropsIndividualFileName(mode, i)    
    _doEnd(mode, i, outFileName, results)

def _louvainClusteringPSpecialWorkFunc(inData):
    mode, i, adjMat = _doStart(inData)
    weighted = False
    g = netgen.makeNetWithLargestLinksAndMST(adjMat, params.pSpecial, weighted=weighted)
    del adjMat #save memory
    gc.collect() #save memory
    results = gencomps.getBestLouvainCommunities(g, nIt=params.louvain_iterations, weighted=weighted)
    outFileName = dataio.getLouvainClusteringIndividualFileName(mode, i)
    _doEnd(mode, i, outFileName, results)