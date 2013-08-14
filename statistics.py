import numpy as np
from numpy.random import RandomState
import multiprocessing

#
# The dataArray often present in this module should be either of shape
# (n1+n2, ) or (n1+n2, nTests)
# 


#SOME FUNCS (AND TRANSFORMATIONS) TO COMPUTE TEST STATISTICS:
##########################################################################################

def meanDifference(dataArray, n1):
    """
    Computes the mean difference (variance normalized mean difference)
    
    Args:
        dataArray:  1 or 2 dimensional array
        n1:         the number of elements in the first group
        
    Returns:
        the mean difference
    """
    return np.average(dataArray[:n1],axis=0)-np.average(dataArray[n1:],axis=0)

def tValue(dataArray, n1):
    """
    Computes the t-value (variance normalized mean difference) for a one or
    two dimensional numpy array
    
    Args:
        dataArray:  array of the values
        n1:         the number of elements in the first group
                    
    Returns:
        the t value
    """
    n2 = dataArray.shape[0]-n1
    var1 = np.var(dataArray[:n1], axis=0) 
    var2 = np.var(dataArray[n1:], axis=0)
    mu1 = np.average(dataArray[:n1],axis=0)
    mu2 = np.average(dataArray[n1:],axis=0)
    return (mu1-mu2)/np.sqrt(var1/n1+var2/n2)

def fisherTrans(corrVals):
    """
    Apply the Fisher transformation to enable stronger effects with high 
    values of correlation. (Variance stabilization) 
    """
    return np.arctanh(corrVals)
    
def simMatrixGroupMeanDiff(matrix, n1):
    return simMatrixGroupMeans(matrix, n1, False)[2]

def simMatrixGroupMeans(matrix, n1, paired):
    """
    Computes the mean of the upper triangle (k=1) for the blocks 
    (0,n-1)*(0,n-1) and (n,2n-1)*(n,2n-1), and their difference (for convenience).
    Also the between groups mean is computed.
    
    If paired setup, the between group diagonal elements are neglected
    i.e. removing the effect of same individuals within control/treatment.
    """
    n2 = matrix.shape[0]-n1
    indices1 = np.triu_indices(n1,k=1)
    indices2base = np.triu_indices(n2, k=1)
    indices2I = indices2base[0].copy()+n1
    indices2J = indices2base[1].copy()+n1
    indices2 = (indices2I, indices2J)
    mean1 = np.average(matrix[indices1])
    mean2 = np.average(matrix[indices2])
    #between group similarities:
    incidenceMatrix12 = np.ones( (n2, n1) )
    if paired:
        incidenceMatrix12 -= np.eye(n1)
    betweenI = incidenceMatrix12.nonzero()[0].copy()
    betweenI += n1
    betweenJ = incidenceMatrix12.nonzero()[1].copy()
    betweenIndices = (betweenI, betweenJ)
    meanBetween = np.average(matrix[betweenIndices])
    return mean1, mean2, mean1-mean2, meanBetween


### INTERNAL FUNCTIONS TO COMPUTE PERMUTATIONS:
##########################################################################################

def _getPermutation(pairedStudy, n1, n2, rng, i):
    """
    Helper function to compute a permutation with all the switches.
    
    Args:
        pairedStudy: True or False
        n1: number of instances in the first group
        n2: number of instances in the second group
        rng: A random number generator (None for deterministic permutations)
        i: The "index" of the permutation, when perm_samples = 'all' is used
    """
    if pairedStudy==True:
        if rng == None: #go through full permutation dist.
            return _getPairedPermutation(n1,i)
        else: # a paired setting but a random permutation
            return _getRandomPairedPermutation(n1+n2, rng)
    else: #not a paired setting
        return rng.permutation(n1+n2)    
    
def _getRandomPairedPermutation(nTot, rng):
    n = nTot/2
    perm = np.arange(0,nTot, dtype = np.int32)
    rands = rng.randint(0,2,n)
    perm[:n] += rands*n
    perm[n:] -= rands*n
    return perm
    
def _getPairedPermutation(n,i):
    permBase = np.zeros(n, dtype = np.int32)
    #compute the bin repr.    
    for j in np.linspace(n-1,0,n).astype(np.int32):
        c = i//2**j
        permBase[j] = c      
        i -= 2**j*c
    perm = np.arange(0,2*n, dtype = np.int32)
    perm[:n] += permBase*n
    perm[n:] -= permBase*n
    return perm 
    
    

### TRADITIONAL PERMUTATION TESTS:
##########################################################################################



def permutationTest(data, statFunc, permuteFunc, n1, n2, pairedStudy, perm_samples, seed=-1):
    """
    Performs a permutation test with user specified test statistics 
    yielding an estimate of the pvalue, and the test statistic.
    
    Args:
        statfunc:       one of the following
                            statistics.meanDifference
                            statistics.tValue
        dataArray:      1d or 2d array of values (treatment/patients+control)
        n1:             the number of subjcets in the first group
                        (13 treatment -> n=13)
        n2:             # of subjects in the 2nd group
        perm_samples:   number of permutations as int 
                        OR 'all' (possible only when pairedStudy = True)
        pairedStudy:    True/False: is the test paired or not?
        seed:           seed value for the random number generator
                        (sometimes needed for parallellism)
                        
    Returns:
        Two-sided pvalues and the test statistic (tvalue/mean difference)
    """
    #all pairs are paired -> change perm_samples
    if perm_samples=='all':
        perm_samples = 2**n1 # (n-1): count everything twice for simplicity (a bit foolish)
        rng = None #no rng needed if deterministic permutations
    else:
        rng = _getRandomState(seed)
        
    origStat = statFunc(data, n1)
    leftNs = np.zeros(np.shape(origStat))
    rightNs = np.zeros(np.shape(origStat))
    for i in range(0,perm_samples):
        perm = _getPermutation(pairedStudy, n1, n2, rng, i)
        permdata = permuteFunc(data, perm) 
        stats = statFunc(permdata,n1)
        leftNs += (stats <= origStat)
        rightNs += (stats >= origStat)    
    tailNs = np.minimum(leftNs, rightNs)
    #multiply by two to get two-sided p-values:
    if rng == None:
        return 2.*tailNs/perm_samples, origStat 
    else:
        pvals = (tailNs+1.0)/(perm_samples+1.0)
        return 2*pvals, origStat
    
def arrayPermute(dataArray, perm):
    return dataArray[perm]
    
def matrixPermute(matrix, perm):
    return matrix[perm,:][:,perm]
    
def meanDifferencePermutationTest(dataArray, n1, n2, pairedStudy, perm_samples, seed=-1, nCpus=1):
    return parallelPermutationTest(dataArray, meanDifference, arrayPermute, n1, n2, pairedStudy, perm_samples, seed, nCpus)
    
def tValuePermutationTest(dataArray, n1, n2, pairedStudy, perm_samples, seed=-1, nCpus=1):
    return parallelPermutationTest(dataArray, tValue, arrayPermute, n1, n2, pairedStudy, perm_samples, seed, nCpus)

def groupMeanDiffSimMatrixPermutationTest(simMatrix, n1, n2, pairedStudy, perm_samples, seed=-1):
    "No parallelism available - nor needed (usually at least)"
    return permutationTest(simMatrix, simMatrixGroupMeanDiff , arrayPermute, n1, n2, pairedStudy, perm_samples, seed)
    
def parallelPermutationTest(dataArray, statFunc, permuteFunc, n1, n2, pairedStudy, perm_samples, seed=-1, nCpus=1):
    """
    Returns the permutation p-values and the values of the test statistic
    """
    #in case of no parallellism
    if nCpus == 1:
        return permutationTest(dataArray, statFunc, permuteFunc, n1, n2, pairedStudy, perm_samples, seed=-1, nCpus=1)
    
    #if parallelism:    
    dataArraySlices = np.array_split(dataArray, nCpus, axis=1)
    #create input arguments
    inputArgsList = []

    for i in range(0,nCpus):
        inputArgsList.append( (dataArraySlices[i], statFunc, permuteFunc, n1, n2, pairedStudy, perm_samples, seed) )
    
    #run the estimation
    pool = multiprocessing.Pool(processes=nCpus)
    result = pool.map_async(_parallelPermutationTestHelper, inputArgsList, chunksize=1)
    #hack (to enable keyboard interruption)
    outputs = result.get(31536000)
    pvals = outputs[0][0]
    origStat = outputs[0][1]
    for i in range(1,len(outputs)):
        pvals = np.hstack( (pvals, outputs[i][0]) )
        origStat = np.hstack( (origStat, outputs[i][1]) )
    return pvals, origStat
    
def _parallelPermutationTestHelper(inputArgs):
    """ This function needs to be outside of the parallelPermutationTestHelper """    
    return permutationTest(*inputArgs)


### (p)FDR ESTIMATION WITH PERMUTATIONS:
##########################################################################################



def permTValueFDREstimate(dataArray, n1, n2, pairedStudy, perm_samples, symThresholds=[], asymThresholds=[], seed=1234, nCpus = 1):
    """
    Perform a sam-like permutation test for the given data using predefined 
    thresholds. For each treshold, the number of rejections for the null
    hypothesis (permutation distribution) and the estimated fdr is returned.
    
    Args:
        dataArray:      an np.array for which to compute the test
        n1:             the number of subjects in the first group
        n2:             the number of subjects in the 2nd grouop        
        pairedStudy:    True/False
        perm_samples:   an int or 'all' (if pairedStudy = True)        
        symThresholds:  list of positive! threshold (t-)values
        asymThresholds: list of positive and negative threshold (t-)values
        seed:           seed for the random number generator
        
    Returns:
        fdrs:   
            an np.array containing the estimated fdr values for each 
            treshold
        origRejs:   
            an np.array containing the number of rejections of the orig. 
            data for each treshold (for convenience)
        avgNullRejs:
            an np.array containing the average number of rejections for 
            coming from the permutation distribution for each tresh.
        thresholds:       
            the thresholds used (compute rejections yourself)
    """
    dataArraySlices = np.array_split(dataArray, nCpus, axis=1)
    #create input arguments
    inputArgsList = []
    for i in range(0,nCpus):
        inputArgsList.append( (dataArraySlices[i], n1, n2, pairedStudy, perm_samples, symThresholds, asymThresholds, seed) )

    #run the estimation
    pool = multiprocessing.Pool(processes=nCpus)
    result = pool.map_async(_permTValueFDREstimateHelper, inputArgsList, chunksize=1)
    #hack (to enable keyboard interruption)
    outputs = result.get(31536000)

    totOrigRejsSym = np.zeros( len(symThresholds), dtype=np.int64)
    totNullRejsSym = np.zeros( len(symThresholds), dtype=np.int64)
    totRejsPositiveSym = np.zeros( (perm_samples, len(symThresholds)),dtype=np.int64)
        
    totNullRejsAsym = np.zeros( len(asymThresholds), dtype=np.int64)
    totOrigRejsAsym = np.zeros( len(asymThresholds), dtype=np.int64)
    totRejsPositiveAsym = np.zeros( (perm_samples, len(asymThresholds)), dtype=np.int64)
    
    if perm_samples == 'all':
        perm_samples = 2**n1
   
    for output in outputs:
        totOrigRejsSym += outputs[0]
        totNullRejsSym += outputs[1]
        totRejsPositiveSym += outputs[2]
        
        totNullRejsAsym += outputs[3]
        totOrigRejsAsym += outputs[4]
        totRejsPositiveAsym += outputs[5]

    #count number of positive indices    
    totRejsPositiveSym = np.sign(totRejsPositiveSym)
    posRejsByThreholdsSym = np.sum(totRejsPositiveSym, axis=0) #used for pFDR
    totRejsPositiveAsym = np.sign(totRejsPositiveAsym)
    posRejsByThreholdsAsym = np.sum(totRejsPositiveAsym, axis=0) #used for pFDR
    
    #compute FDRs
    FDRsSym = totNullRejsSym/(perm_samples*totOrigRejsSym)
    FDRsSym[np.isposinf(FDRsSym)] = 0 #change to infs to zero (FDR definition)
    FDRsAsym = totNullRejsAsym/(perm_samples*totOrigRejsAsym)
    FDRsAsym[np.isposinf(FDRsAsym)] = 0 #change to infs to zero (FDR definition)
    
    #compute pFDRs
    #pFDR not defined when totOrigRejs == 0 -> inf
    pFDRsSym = totNullRejsSym/(posRejsByThreholdsSym*totOrigRejsSym)
    pFDRsAsym = totNullRejsSym/(posRejsByThreholdsAsym*totOrigRejsAsym)
    
    return FDRsSym, pFDRsSym, FDRsAsym, pFDRsAsym

def _permTValueFDREstimateHelper(inputArgs):
    """ This function needs to be outside of the parallelPermFDREstimate """
    return _permTValueFDREstimateSlave(*inputArgs)


def _permTValueFDREstimateSlave(dataArray, n1, n2, pairedStudy, perm_samples, symThresholds=[], asymThresholds=[], seed=-1):
    """
    See permTValueFDREstimate for actual use.
    """
    #all pairs are paired -> change perm_samples
    if perm_samples=='all':
        perm_samples = 2**n1 # vs. (n1-1): count everything twice for simplicity (a bit foolish..)
        rng = None #no rng needed if deterministic permutations
    else:
        rng = _getRandomState(seed)
        
    origRejsSym = np.zeros( len(symThresholds) , dtype=np.int64)
    nullRejsSym = np.zeros( len(symThresholds) , dtype=np.int64)
    rejsPositiveSym = np.zeros( (perm_samples, len(symThresholds)), dtype=np.int64 )
    
    
    nullRejsAsym = np.zeros( len(asymThresholds), dtype=np.int64)
    origRejsAsym = np.zeros( len(asymThresholds), dtype=np.int64)
    rejsPositiveAsym = np.zeros( (perm_samples, len(asymThresholds)), dtype=np.int64)
    
    asymThresholdsSigns = np.sign(asymThresholds)

    #compute original rejections    
    #symmetric (two-sided) thresholds    
    origStats = tValue(dataArray)
    for i, thresh in enumerate(symThresholds):
        origRejsSym[i] = np.sum(np.abs(origStats) >= thresh)
    #asymmetric (one-sided) thresholds:
    for i, thresh in enumerate(asymThresholds):
        sign = asymThresholdsSigns[i]
        origRejsAsym[i] = np.sum(sign*origStats >= sign*thresh)
        
    #do the permutations
    for i in range(0,perm_samples):
        perm = _getPermutation(pairedStudy, n1, n2, rng, i)
        permdata = dataArray[perm]
        stats = tValue(permdata,n1)
        #symmetric (two-sided) thresholds:
        for j, thresh in enumerate(symThresholds):
            rejs = np.sum(np.abs(stats) >= thresh)
            nullRejsSym[j] += rejs
            rejsPositiveSym[i][j] = np.sign(rejs) #for pFDR
        #asymmetric (one-sided) thresholds:
        for j, thresh in enumerate(asymThresholds):
            sign = asymThresholdsSigns[i]
            rejs = np.sum(sign*stats >= sign*thresh) 
            nullRejsAsym[j] += rejs    
            rejsPositiveAsym[i][j] = np.sign(rejs) #for pFDR
    return origRejsSym, nullRejsSym, rejsPositiveSym, origRejsAsym, nullRejsAsym, rejsPositiveAsym


### P-VALUE CORRECTIONS ETC.:
##########################################################################################

def getFdrTresholdBH(q, pvalues):
    """
    Computes the fdr treshold based on the given p-values and the q-value.
    
    Args:
        q: the estimated rate of false discoveries
        pvalues: The (non-adjusted) p-values obtained from the individual 
                 tests.
                 type preferably flat numpy array.
    
    Returns:
        A single floating point number describing the treshold p-value.
        If all observations should be accepted, returns 1.
        
    References:
        Wikipedia, False discovery rate:
        http://en.wikipedia.org/wiki/False_discovery_rate
        
    Notes:
        Assumes independence for the individual p-values.
        Fast with even 1e6 pvalues (or more...)
    """
    pvalscopy = pvalues.flatten()   #pvals as a vector, flatten = copy
    pvalscopy.sort()                #sorted from smallest to largest
    m = len(pvalscopy)              #number of pvals
    compvec = np.arange(1,m+1)*(q/float(m)) #to which we compare
    fulfilling_indices = np.nonzero(pvalscopy <= compvec )
    if len(fulfilling_indices[0]) == 0:
        return 0.0
    return compvec[np.max(fulfilling_indices)]


    
def pValuesToPFDRWithMatlab(pValues):
    """ 
    Simple matlab wrapper for mafdr (pFDR+estimation of null by Storey)

    Params:
        pValues: a numpy (1d) array containing pValues for which pFDR and/or
                    q-values should be estimated
                    
    Returns:
        fdrs: a computed pFDR level for each p-value
        qvals: a q-value for each p-value (~prob of having null if rejection 
                limit was determined by the p-value)
        pi0: the estimated proportion of the null distribution 
                (when two-dist. model is assumed)                    
    """
    from mlabwrap import mlab #start matlab only if needed
    fdrs, qvals, pi0 = mlab.mafdr(pValues, nout=3)
    fdrs = fdrs.ravel()
    qvals = qvals.ravel()
    pi0 = float(pi0)
    return fdrs, qvals, pi0


### BOOTSTRAP COMPUTATIONS
##########################################################################################

def bootstrapMeanConfInterval(data, nIt, percent):
    """
    Compute the mean confidence interval using bootstrap resampling
    
    Args:
        data:   the (1d) 2d array for which to compute the mean confidence interval
        nIt:    the number of resamplings to 
        percent: the total coverage in percentages (symmetric)
    Returns:
        The percentiles of the data [prop, percentile]
    """
    tail = (100-percent)/2.
    percentiles = [tail, 100-tail]
    indices = np.random.randint(0,data.shape[-1], (nIt, data.shape[-1]) )
    avgs = np.average(np.array(data)[...,indices], axis=-1)
    return np.array(np.percentile(avgs, percentiles, axis=-1)).T
    
def bootstrapMeanConfIntervalForSimMatrix(matrix, n1):
    """
    
    """
    n2 = matrix.shape[0]-n1
    indices1 = np.triu_indices(n1,k=1)
    indices2base = np.triu_indices(n2, k=1)
    indices2I = indices2base[0].copy()+n1
    indices2J = indices2base[1].copy()+n1
    indices2 = (indices2I, indices2J)
    values1 = matrix[indices1]
    values2 = matrix[indices2]
    mean1 = np.average(values1)
    mean2 = np.average(values2)
    return np.array( [mean1, mean2] ), bootstrapMeanConfInterval(np.vstack( (values1, values2) ) )
    
### OTHER:
##########################################################################################  


def _getRandomState(seed=None):
    """
    Get a numpy.random.RandomState object with the given seed.
    If no seed is given or seed = None, a RandomState object with default 
    numpy seeding method (from time) is returned.
    
    Args:
        seed: an int or None
        
    Returns:
        The RandomState object (a random number generator)
    """
    if seed is not None:
        return RandomState(seed)
    else:
        return RandomState()
        
def tValueThresholds(dataArray, n, percents, asym=True):
    """
    asym==True:
        Get left and right tail tvalue thresholds for given percentiles
    asym==False:
        Get tValue treholds for varying percentage of tail links.
        E.g. percents = [5] returns a t-value treshold which roughly counts for 
        taking 2*5 = 10 percent of the links (two tails)
    """
    tVals = tValue(dataArray, n)
    lower = np.percentile(tVals, list(percents) ) 
    upper = np.percentile(tVals, list(100-np.array(percents))) 
    if asym:
        return np.hstack( (lower, upper) )
    else:
        return np.average(np.vstack( (upper,np.abs(lower)) ),0)


### OLD/NOT IN USE:
##########################################################################################

#        
#def _checkPositiveness(thresholds):
#    for threshold in thresholds:
#        if threshold < 0:
#            raise StandardError("A threshold value encountered was less than zero")

#def _getNTests(dataArray):
#    """
#    Resolves the number of tests to be performed based on the given np.array
#    
#    Args:
#        dataArray:  the np.array for which the number of tests should be 
#                    computed
#                    
#    Raises:
#        StandardError if the matrix dimensionality is greater than 2.
#        
#    Returns:
#        The length of the second dimension of the array if input is 2d or 
#        1 if 1-dimensional array.
#    
#    """
#    #i.e. take care of the dimensions of the given array
#    nDim = dataArray.ndim
#    if nDim > 2:
#        raise StandardError("Invalid sized matrix sent to statistical testing")
#    elif nDim == 2:
#        nTests = dataArray.shape[1]
#    else:
#        nTests = 1
#    return nTests
#        

    