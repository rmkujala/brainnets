import numpy as np
import scipy.io
import igraph
import itertools
from collections import defaultdict

def getMutualInformationCommunities(dataType='Best'):
    '''
    Get the similarity between modules across the subjects
    dataType: 'Best', 'Consensus'
    '''
    if dataType=='Best':
        dirName='BestModules/'
        filePrefix='bestModule_subject_'
        outFileName=dirName+'bestModuleComparision.mat'
    elif dataType=='Consensus':
        dirName='ConsensusModules/'
        filePrefix='consensusModule_subject_'
        outFileName=dirName+'consensusModuleComparision.mat'

    nSub=26 # Number of subjects
    nPercentage=50 # Number of % where the properties are calculated
    metrices=['vi','nmi','split-join','rand','adjusted_rand']
    mArray={}
    for m in metrices:
        mArray[m]=np.zeros((nSub,nSub,nPercentage))

    comms=defaultdict(list)
    for s in range(nSub):
        commInFileName=dirName+filePrefix+str(s+1)+'.txt'
        with open(commInFileName) as fp:
            for line in fp:
                comms[s].append(map(int,line.split()))

    for p in range(nPercentage):
        print p
        for m in metrices :
            for s1,s2 in itertools.combinations(range(nSub),2):
                mArray[m][s1,s2,p]=igraph.compare_communities(comms[s1][p],comms[s2][p],method=m)
                mArray[m][s2,s1,p]=mArray[m][s1,s2,p]

    mArray['split_join']=mArray.pop('split-join')
    scipy.io.savemat(outFileName,mArray)

def getMutualInformationBestModules():
    '''
    Get the similarities between modules across the subjects 
    - For Best modules
    '''
    getMutualInformationCommunities('Best')
    
def getMutualInformationConsensusModules():
    '''
    Get the similarities between modules across the subjects 
    - For Best modules
    '''
    getMutualInformationCommunities('Consensus')

def getNodesPearsonPhi(dataType='Best',pLinks=2):
    '''
    Get the similarity between modules across the subjects
    dataType: 'Best', 'Consensus'
    pLinks = % of Links considered
    '''
    print pLinks
    if dataType=='Best':
        dirName='BestModules/'
        filePrefix='bestModule_subject_'
        outFileName=dirName+'bestModule_pearsonPhi_pLinks_'+str(pLinks)+'.mat'
    elif dataType=='Consensus':
        dirName='ConsensusModules/'
        filePrefix='consensusModule_subject_'
        outFileName=dirName+'consensusModule_pearsonPhi_pLinks_'+str(pLinks)+'.mat'

    nSub=26 # Number of subjects
    comms=[]
    for s in range(nSub):
        commInFileName=dirName+filePrefix+str(s+1)+'.txt'
        c=np.loadtxt(commInFileName,skiprows=pLinks-1,dtype=int)[0,:]
        comms.append(c)

    nNodes=comms[0].shape[0]
    mArray=np.zeros((nSub,nSub,nNodes))

    for node in range(nNodes):
        if node%100==0:
            print node
        for s1,s2 in itertools.combinations(range(nSub),2):
            mArray[s1,s2,node]=np.corrcoef(comms[s1]==comms[s1][node],comms[s2]==comms[s2][node])[0][1]
            mArray[s2,s1,node]=mArray[s1,s2,node]

    scipy.io.savemat(outFileName,{'phi':mArray})

def getNodesPearsonPhiBestModules():
    '''
    Get the pearson phi across the subjects 
    - For Best modules
    '''
    for p in [2,5,10]:
        getNodesPearsonPhi(dataType='Best',pLinks=p)
    
def getNodesPearsonPhiConsensusModules():
    '''
    Get the pearson phi across the subjects 
    - For Consensus modules
    '''
    for p in [2,5,10]:
        getNodesPearsonPhi(dataType='Consensus',pLinks=p)

def testSignificantDifference(nPermutations=1000):

    mat=scipy.io.loadmat('BestModules/bestModuleComparision.mat')
    nmi=mat['nmi'][:,:,1] # nmi at 2%
    #nmi=mat['adjusted_rand'][:,:,1] # nmi at 2%

    randNMI=np.zeros(nPermutations)
    for nPerm in range(nPermutations):
        idx=np.random.permutation(26)
        randNMI[nPerm]=nmi[idx,:][:,idx][:13,13:].mean()

    meanNMI=nmi[:13,13:].mean()
    print np.sum(randNMI<meanNMI)

def testSignificantDifference2(nPermutations=1000):

    mat=scipy.io.loadmat('BestModules/bestModuleComparision.mat')
    nmi=mat['nmi'][:,:,1] # nmi at 2%
    #nmi=mat['adjusted_rand'][:,:,1] # nmi at 2%

    randNMI=np.zeros(nPermutations)
    diagIdx=np.triu_indices(13,1)
    for nPerm in range(nPermutations):
        idx=np.random.permutation(26)
        randNMI[nPerm]=nmi[idx,:][:,idx][:13,13:].mean()-nmi[idx,:][:,idx][diagIdx].mean()

    meanNMI=nmi[:13,13:].mean()-nmi[diagIdx].mean()
    print np.sum(randNMI<meanNMI)
