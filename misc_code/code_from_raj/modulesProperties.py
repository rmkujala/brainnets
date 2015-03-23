#!/usr/bin/env python
# Various Properties of modules, subgraph or network
# e.g., modularity, coherence, intensity 
import numpy as np
import time
import itertools

def getIntensityCoherenceSubgraph(weights):
    '''
    Get the intensity and coherence of the subgraph
    PHYSICAL REVIEW E 71, 065103R 2005
    Input:
    List of all the weights in the subgraph
    '''
    sWeight = sum(weights)
    pWeight = reduce(mul,weights)
    k=len(weights)
    intensity=pWeight**(1./k)
    coherence=(intensity*k)/sWeight
    return intensity, coherence

def getIntensityCoherence(commProps,nModules):
    '''
    Get the intensity and coherence 
    PHYSICAL REVIEW E 71, 065103R 2005
    '''
    intensity=np.zeros(nModules)
    coherence=np.zeros(nModules)
    idx=0
    for k, sWeight, pWeight, in commProps:
        intensity[idx]=pWeight**(1./k)
        coherence[idx]=(intensity[idx]*k)/sWeight
        idx += 1
    return intensity, coherence

def getConsensusMatrix(nodeList,cMatrix):
    '''
    Get the consensus matrix
    '''
    x=np.array(nodeList)
    for comm in range(x.max()+1):
        idx=(x==comm).nonzero()[0]
        for i, j in itertools.combinations(idx,2):
            cMatrix[i,j]+=1
    return cMatrix

def getModularCompactness(nodeList,distArray):
    '''
    Average distance (Euclidean) between all possible pair of nodes
    '''
    x=np.array(nodeList)
    nModules=x.max()+1
    commCompactness=np.zeros(nModules, dtype=[('meanEuclDist', 'f4')])
    for comm in range(nModules):
        idx=(x==comm).nonzero()[0]
        commCompactness[comm]['meanEuclDist']=distArray[idx,:][:,idx].mean()
    return commCompactness

def getAverageModularCompactness(commCompactness,commProps):
    '''
    Get the average compactness (weighted by size of the communities)
    '''
    return (commCompactness['meanEuclDist']*commProps['size']).sum()/commProps['size'].sum()

def getThreholdConsensusMatrix(commProps):
    '''
    Threshold - Null modes assumes that the nodes are distributed uniformly
    across modules 
    '''
    moduleSizes=commProps['size']
    nNodes=moduleSizes.sum()
    return (moduleSizes*(moduleSizes-1)).sum()/(nNodes*(nNodes-1.0))

def communityProperties(net, nodeList):
    '''
    (*) Then number of internal links (l = intra-community Links)
    (*) Sum of the weights of the internal links (w)
    (*) The total degree of all nodes in the community (d) 
    (*) The strength of all nodes in the community (s)
    (*) Size of each community
    The number of inter-community Links (= d-l)
    The sum of the weights of the inter-community links (= s-w)
    '''
    # Number of communities
    comms=list(set(nodeList))  
    nComms=len(comms)
    commProps=np.zeros((nComms,), dtype=[('nLinks', 'i4'), ('linkWeights', 'f4'), ('degree','i4'),('strength','f4'),('size','i4')])
  
    # Iterate through each edge and if both the nodes of the edge are in same community increase the intra-links for that community
    for i, j, weight in net.edges:
        if nodeList[i]==nodeList[j]:
          commProps[nodeList[i]]['nLinks'] += 1;
          commProps[nodeList[i]]['linkWeights'] += weight;

    # Determine the sum of the degree of nodes in each community
    for i in net:
        commProps[nodeList[i]]['degree'] += net[i].deg();
        commProps[nodeList[i]]['strength'] += net[i].strength();
        commProps[nodeList[i]]['size'] += 1;
    return commProps

def getModularity(commProps, nLinks, gamma=1): 
    '''
    Modularity of the undirected network  
    Input:
    For each module a tuple with inside Links and total degree 
    nLinks : Total links
    Optional:
    gamma : Resolution parameter(Float, default=1)
    '''
    modularity=0
    for links, deg in commProps[['nLinks','degree']]:
        modularity += links/float(nLinks) - gamma*((deg/(2.0*nLinks))**2)
    return modularity

def getInLinksAvgWeight(commProps):
    '''
    Average weight of links inside community
    '''
    return commProps['linkWeights'].sum()/commProps['nLinks'].sum()

def getDensityInLinks(commProps):
    '''
    Inter-modular link density (weighted by size)
    \sum n_i/N l_i/(n_i*(n_i-1))
    '''
    idx=(commProps['size']!=1)
    return ((2.0*commProps['nLinks'][idx])/(commProps['size'][idx]-1.)).sum()/(commProps['size'][idx]).sum()

def getTotalInLinks(commProps):
    '''
    Number of Inter-modular links
    '''
    return commProps['nLinks'].sum()

def netWeightShuffle(net,weightList=None):
    '''
    Shuffle the weights of the network
    '''
    if weightList==None :
        weightList=np.array(list(net.weights))

    nEdges=weightList.shape[0]
    weightList=weightList[np.random.permutation(nEdges)]
    idx=0
    for node1, node2, weight in net.edges:
        net[node1,node2]=weightList[idx]
        idx += 1
    return net

