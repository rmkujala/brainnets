import scipy.io
import numpy as np
import cPickle as pickle
import os
import brainnets_params as params
	
def loadPickle(fName):
    f = open(fName, "rb")
    data = pickle.load(f)
    f.close()
    return data
    
def savePickle(fName, dataDict, overwrite=True):
    """
    Saves the dataDict to a pickle file.
    If the pickled file exists (with data), the results will just be appended.
    (In order to not lose any data)
    """
    if not overwrite and os.path.exists(fName):
        #Be sure that dataDict is really a dictionary.
        try:
            f = open(fName, "rb")
            data = pickle.load(f)
            for key in dataDict:
                data[key] = dataDict[key] #update old ones
            f.close()
            f = open(fName, "wb") #overwrite the old file
            pickle.dump(data, f)
            f.close()
        except:
            print "Something went wrong while saving/appending pickle data"
            raise
    else: #default behaviour
        f = open(fName, "wb")
        pickle.dump(dataDict, f)
        f.close()

def loadAdjMatrixFromMat(filename, varName = 'adj'):
    """
    Load a correlation matrix using .mat file format 
    """
    assert filename[:-4] == ".mat", "Trying to load incorrect file format (.mat)."
    adjMatrix = scipy.io.loadmat(filename)[varName]
    adjMatrix = np.triu(adjMatrix,1)
    return adjMatrix
        
def texOutputArray(array, outFileName, prec=2):
    """
    Outputs an numpy array to a .tex format that can be inported to a master .tex document. 
    """    
    string = " \\\\\n".join([" & ".join(map( ('{0:.'+str(prec)+'f}').format, line)) for line in array])
    f = open(outFileName, "w")
    f.write(string)
    f.close()


def outPutNodePropToNii(nodeVals, outFileName, roiInfoMFileName, blacklistFileName):
    """
    Calls a matlab function to output the given values of a node property to
    a nii volume.

    Arguments:
        nodeVals: list of values with the original indexing (non blacklisted)
        outFileName: the name of the file into which the volume info is written
                        (should end with a .nii)
        roiInfoFileName: the filename containing the ROI info 
                        (passed directly to matlab)
        blacklist: the list of blacklist nodes
    """
    from mlabwrap import mlab 
    mlab.addpath("/proj/networks/rmkujala/brain_fmri/code/src")
    mlab.data2niivol(nodeVals, blacklistFileName, roiInfoMFileName, outFileName)
    return


def getBlacklistedIndicesFromMat():
    """Returns the blacklist using python indexing (starting from 0)"""
    fName = params.blacklistNodesFileName
    blacklist = scipy.io.loadmat(fName, squeeze_me = True)["blacklist"]
    #in matlab indexing starts from one:
    return np.array(blacklist)-1 

def getOkIndices():
    """Returns the valid nodes, based on the blacklist"""
    nOld=params.mrNVoxels
    blacklist=getBlacklistedIndicesFromMat()
    old_indices = range(0,nOld)
    newToOldMap = [index for index in old_indices if index not in blacklist]
    return newToOldMap

def filterBlacklisted(matrix):
    """Filters a blacklisted matrix"""
    nonBlNodes = getOkIndices()
    filteredmat = matrix[nonBlNodes, :]
    filteredmat = filteredmat[:,nonBlNodes]
    return filteredmat

def expand1DNodeValArrayToOldIndices(values):
    """
    Expands a 1d node val array to indices used before removing the 
    blacklisted nodes
    """
    nOld=params.nVoxels
    newVals = np.zeros(nOld)    
    newVals[getOkIndices()] = values
    return newVals

def loadFilteredAdjacencyMatrix(mode, subject_number):
    """
    Loads blacklist-filtered adjacency matrix into memory and returns it
    
    Args:
        mode: params.mode1 or params.mode2
        subject_number: e.g 1
    Returns:
        The blacklist-filtered matrices in a list
    """
    fName = getCorrMatFname(mode, str(subject_number) )
    return filterBlacklisted(loadAdjMatrixFromMat(fName))
    
def loadFlattenedAdjacencyMatrix(mode, i):
    matrix = loadFilteredAdjacencyMatrix(mode, i)
    twoDindices = np.triu_indices_from(matrix,1)
    return matrix[twoDindices].flatten()
    
def loadFlattenedAdjacencyMatrices():
    """
    Load all mode1 and mode2 correlation matrices.
    """
    #computing these beforehand + pickling -> error with unpickling 
    #(a bug in the pickle library it seems, memory stuff..)
    corrMatrices = []
    print "starting to load correlation matrices..."
    for i in range(1,params.n1+1):
        corrMatrices.append(loadFlattenedAdjacencyMatrix(params.mode1, i))
    for i in range(1,params.n2+1):
        corrMatrices.append(loadFlattenedAdjacencyMatrix(params.mode2, i))
#    twoDindices = np.triu_indices_from(matrices[0], 1)
#    data = []
#    for matrix in matrices:
#        data.append(matrix[twoDindices].flatten())
#    data = np.array(data)
    return np.array(corrMatrices)

def mergeAndLoadNodeProperties():
    """
    Load the individual node properties and combine them to a joint 
    dictionary with structure:
    data[prop][subjectNmovie->rest][nodeIndex] = value
    """
    return mergeAndLoadProperties(loadIndNodeProperties, params.nodeProps)
    
def mergeAndLoadGlobalUWProperties():
    """
    Load the individual global UW properties and combine them to a joint 
    dictionary with structure:
    data[prop][subjectNmovie->rest][p/nodeIndex/linkIndex] = value
    """
    return  mergeAndLoadProperties(loadIndGlobUWProperties, params.globUWProps)
    
def mergeAndLoadGlobalWProperties():
    """
    Load the individual global UW properties and combine them to a joint 
    dictionary with structure:
    data[prop][subjectNmovie->rest][p/nodeIndex/linkIndex] = value
    """
    return mergeAndLoadProperties(loadIndGlobWProperties, params.globWProps)
    
def mergeAndLoadLouvainProperties():
    return mergeAndLoadProperties(loadIndLouvainProperties, params.louvainProps)
    
def mergeAndLoadProperties(indPropLoadFunc, props):
    """
    Load the individual properties (glob UW, glob W, node, link) 
    and combine them to a joint file with structure:
    data[prop][subjectN(movie...rest)][nodeIndex] = value
    """
    dataDict = {}
    counter = 0
    #loop over modes
    for mode in params.modes:
        #loop over subjects        
        for i in range(1,params.mrNsubjs+1):
            indData = indPropLoadFunc(mode, i)
            #loop over properties
            for prop in set(indData.keys()).intersection(props):
                #initialize new properties:                    
                if counter == 0:
                    #l = number of nodes / percetages / or zero (modularity stuff)
                    try:
                        l = len(indData[prop])
                    except:
                        l=1
                    dataDict[prop] = np.zeros( (params.n1+params.n2, l))
                #add the data to the dict:
                dataDict[prop][counter,:] = indData[prop]
            counter += 1
            try:
                dataDict[params.percentages_abbr] = indData[params.percentages_abbr]
            except:
                pass
    return dataDict


def loadIndNodeProperties(mode, subjectNumber):
    return loadPickle(params.getNodePropsIndividualFileName(mode, subjectNumber) )

def loadIndGlobUWProperties(mode, subjectNumber):
    return loadPickle(params.getGlobalUWPropsIndividualFileName(mode, subjectNumber) )

def loadIndGlobWProperties(mode, subjectNumber):
    return loadPickle(params.getGlobalWPropsIndividualFileName(mode, subjectNumber) )
    
def loadIndLouvainProperties(mode, subjectNumber):
    return loadPickle(params.hetLouvainClusteringIndividualFileName(mode, subjectNumber) )
    
###
### Filenaming conventions:
###

def getCorrMatFname(mode, i):
    """
    Returns the filename of the correlation matrix    
    """
    return params.inputDir+mode+"_"+str(i)+".mat"

def getGlobalUWPropsIndividualFileName(mode, subjectNumber):
    return params.outputDir + "globalUWProps_"+indFileEnding(mode, subjectNumber)
        
def getGlobalWPropsIndividualFileName(mode, subjectNumber):
    return params.outputDir + "globalWProps_"+indFileEnding(mode, subjectNumber)
    
def getNodePropsIndividualFileName(mode, subjectNumber):
    return params.outputDir + "nodeProps_"+indFileEnding(mode, subjectNumber)
        
def getLinkPropsIndividualFileName(mode, subjectNumber):
    return params.outputDir + "linkProps_"+indFileEnding(mode, subjectNumber)

def getLouvainClusteringIndividualFileName(mode, subjectNumber):
    return params.outputDir + "louvainClustering_"+indFileEnding(mode, subjectNumber)

def indFileEnding(mode, subjectNumber):
    return mode+"_%d.pkl" %(subjectNumber)
    
    
def getGlobalUWStatsFilename():
    return params.outputDir+ "globalUWStats.pkl"

def getGlobalWStatsFilename():
    return params.outputDir + "globalWStats.pkl"
    

#p-values    
def getNodePValStatsFilename():
    return params.outputDir + "nodePVals.pkl"
    
def getCorrPValStatsFilename():
    return params.outputDir + "corrPVals.pkl"

    
#dependent correlation stats
def getNodeFDRStatsFilename(asym=False):
    baseName = "nodeFDRStats"
    if asym:
        baseName+="_asym"        
    return params.outputDir + baseName+".pkl"
    
def getLinkFDRStatsFilename(asym=False):
    baseName = "linkFDRStats"
    if asym:
        baseName+="_asym"        
    return params.outputDir + baseName+".pkl"
    
def getCorrFDRStatsFilename(asym=False):
    baseName = "corrFDRStats"
    if asym:
        baseName+="_asym"        
    return params.outputDir + baseName+".pkl"

    
