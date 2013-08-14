### This module is independent of 
### all other modules but settings and netgen
### (no params settings directly used here)



import settings
import igraph
import os
import numpy as np
import subprocess
import sys
import netgen


def getGlobalUWProps(net, props=None):
    """
    Compute the global _unweighted_ properties determined by the argument 
    props.
    
    Args:
        net:    an igraph.Graph instance
        props:  a list of global unweighted properties as shown in fmrifmrisettings.py
                e.g. ["acl", "gcl", "ass"]
                If props are not given, all globalUWprops listed in 
                fmrisettings.py are computed.
                
    Returns:
        measures: a list, where the elements are in the order of the properties
    """
    if props == None:
        props = settings.globUWProps
        
    resultDict = {} #add results in a dict, unwrap later
    if "acl" in props:
        resultDict["acl"] = net.transitivity_avglocal_undirected(mode='zero')
    if "gcl" in props:
        resultDict["gcl"] = net.transitivity_undirected()
    if "apl" in props:
        resultDict["apl"] = net.average_path_length()
    if "ass" in props:
        resultDict["ass"] = net.assortativity_degree()
    if "maxk" in props:
        resultDict["maxk"] =  max(net.degree(net.vs))
    if "mksh" in props:
        resultDict["mksh"] = np.max(net.shell_index())
    return _unwrap(resultDict, props)
    
def getGlobalWProps(net, props=None):
    """
    Compute the global _weighted_ (=correlation) properties determined by the 
    argument props.
    
    Args:
        net:    a weighted igraph.Graph instance (weight=correlation)
        props:  a list of global weighted properties as shown in fmrisettings.py
                e.g. ["wpl", "wass"]
                If props are not given, all globalWprops listed in 
                fmrisettings.py are computed.
                
    Returns:
        resultdict: a dictionary, where key is the property (e.g. "wpl") and
                    the value is the computed measure
    """
    if props == None:
        props = settings.globWProps
        
    resultDict = {} #add results in a dict, unwrap later
    if "wapl" in props:
        swpls= np.array(net.shortest_paths(weights=1./np.array(net.es["weight"])))
        #shortest weighted path lengths
        resultDict["wapl"] = np.average(swpls[np.triu_indices_from(swpls,1)])
    if "maxs" in props:
        resultDict["maxs"] = np.max(net.strength(net.vs, weights=net.es["weight"]))
    if "wcl" in props:
        resultDict["wcl"] = net.transitivity_avglocal_undirected(mode="zero")
    if "wass" in props:
        resultDict["wass"] = net.assortativity(net.strength(weights=net.es['weight']), directed=False)
    return _unwrap(resultDict, props)
    
def getGlobalPropsForPRange(corrMat, pRange, props, weighted):
    resultDict = {}
    for prop in props:
        resultDict[prop] = np.zeros( len(pRange) )
    for i, p in enumerate(pRange):
        print p
        net = netgen.makeNetWithLargestLinksAndMST(corrMat, p, weighted=False)
        if weighted == True:
            results = getGlobalWProps(net, props)
        if weighted == False:
            results = getGlobalUWProps(net, props)
        for j, prop in enumerate(props):
            resultDict[prop][i] = results[j]
    resultDict[settings.percentages_abbr] = pRange
    return resultDict

#def getGlobalUWPropsForPRange(corrMat, pRange, props):
#    resultDict = {}
#    for prop in props:
#        resultDict[prop] = np.zeros( len(pRange) )
#    for i, p in enumerate(pRange):
#        print p
#        net = netgen.makeNetWithLargestLinksAndMST(corrMat, p, weighted=False)
#        results = getGlobalUWProps(net, props)
#        for j, prop in enumerate(props):
#            resultDict[prop][i] = results[j]
#    resultDict[settings.percentages_abbr] = pRange
#    return resultDict
#    
#def getGlobalWPropsForPRange(corrMat, pRange, props):
#    resultDict = {}
#    for prop in props:
#        resultDict[prop] = np.zeros( len(pRange) )
#    for i, p in enumerate(pRange):
#        net = netgen.makeNetWithLargestLinksAndMST(corrMat, p, weighted=True)
#        results = getGlobalWProps(net, props)
#        for j, prop in enumerate(props):
#            resultDict[prop][i] = results[j]
#    resultDict[settings.percentages_abbr] = pRange
#    return resultDict

    
    
def getNodeProps(net, props=None):
    """
    Compute the node level properties determined by the argument 
    props.
    
    Args:
        net:    a weighted igraph.Graph instance
        props:  a list of node weighted properties as shown in fmrisettings.py
                e.g. ["s", "k"]
                If props are not given, all nodeProps listed in 
                fmrisettings.py are computed.
                
    Returns:
        resultdict: a dictionary, where key is the property (e.g. "s") and
                    the value is a np.array containing the value of the 
                    measure for each node.
    """
    if props == None:
        props = settings.nodeProps
    
    resultDict = {} #add results in a dict, unwrap later
    
    if "k" in props:
        resultDict["k"] = net.degree(net.vs)
    if "s" in props:
        resultDict["s"] = net.strength(weights=net.es["weight"])
    if "bc" in props:
        resultDict["bc"] = net.betweenness()
    if "wbc" in props:
        resultDict["wbc"] = net.betweenness(weights=net.es["weight"])
    if "ksh" in props:
        resultDict["ksh"] = net.shell_index()
        
    return _unwrap(resultDict, props)


def getLinkProps(net, props=None):
    """
    Compute the node level properties determined by the argument 
    props.
    
    Args:
        net:    a weighted igraph.Graph instance
        props:  a list of node weighted properties as shown in fmrisettings.py
                e.g. ["bc"]
                If props are not given, all linkProps listed in 
                fmrisettings.py are computed.
                
    Returns:
        resultdict: a dictionary, where key is the property (e.g. "bc") and
                    the value is a np.array containing the value of the 
                    measure for each node.
    """
    if props == None:
        props = settings.linkProps
    resultDict = {} #add results in a dict, unwrap later
    
    if "ebc" in props:
        resultDict["ebc"] = net.edge_betweenness()
    if "webc" in props:
        resultDict["webc"] = net.edge_betweenness(weights=net.es["weight"])
    
    return _unwrap(resultDict, props)


def _unwrap(resultDict, props):
    """
    Unwraps the results a dict to a list
    """
    resultList = []
    for prop in props:
        resultList.append(resultDict[prop])
    return resultList
    

def getNodePropsForP(corrMat, p, props):
    """ Get node properties for a specific value of p """
    resultDict = {}
    net = netgen.makeNetWithLargestLinksAndMST(corrMat, p, weighted=True)
    results = getNodeProps(net, props)
    for i, prop in enumerate(props):
        resultDict[prop] = results[i]
    resultDict[settings.percentages_abbr] = p
    return resultDict
    
def getLinkPropsForP(corrMat, p, props):
    resultDict = {}
    net = netgen.makeNetWithLargestLinksAndMST(corrMat, p, weighted=True)
    results = getLinkProps(net, props)
    for i, prop in enumerate(props):
        resultDict[prop] = results[i]
    resultDict[settings.percentages_abbr] = p
    return resultDict
    
    
    
    
    
    
    
### Module level stuff:
##############################################################################

def getLouvainCommunities(graph, weighted=False):
    """
    Obtain clusters using the Louvain code by the original authors + Raj's randomization.
    (the igraph code seems to give deterministic (=same) results for each iteration)
    
    Author: Raj, adapted by Rainer
    """
    baseName = "/tmp/louvain" + str(os.getpid()) #enables parallelism

    #write graph as edge list
    with open(baseName+".edg", "w") as f:
        if weighted:
            for e in graph.es:
                strToWrite = str(e.source) + " " + str(e.target) + " " + str(e["weight"]) + "\n"
                f.write(strToWrite)
        else:
            for e in graph.es:
                strToWrite = str(e.source) + " " + str(e.target) + "\n"
                f.write(strToWrite)
        
    #convert to bin format
    if weighted == True:
        convertRunstr=settings.louvainConvert+' -i '+baseName+".edg"+' -o '+baseName+'.bin' + ' -w ' + baseName+".weights"
    else:
        convertRunstr=settings.louvainConvert+' -i '+baseName+".edg"+' -o '+baseName+'.bin'
    subprocess.call(convertRunstr.split())
    
    # Detect the communities
    runstr = settings.louvainCommunity + " " +baseName+'.bin'+' -l -1 -r 1'
    if weighted==True :
        runstr += ' -w '+baseName+".weights"
    with open(baseName+'.tree','w') as fp:
        p = subprocess.Popen(runstr.split(),stdout=fp,stderr=subprocess.PIPE)
    _,q=p.communicate() #modularity score

    # Get the number of hierarchical level
    runstr = settings.louvainHierarchy + " " + baseName+".tree"
    p = subprocess.Popen(runstr.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    hLevels,_=p.communicate()
    hLevel=int(hLevels.split('\n')[0].split(':')[1])-1
    # Nodes in the final hierarchical level
    runstr = settings.louvainHierarchy + " " + baseName+".tree" + ' -l ' + str(hLevel) 
    p = subprocess.Popen(runstr.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    nodeOutput,_=p.communicate()
    
    # Read the community
    nodeList=[]
    for line in nodeOutput.split('\n')[:-1]: 
        nodeList.append(int(line.split()[1]))
    
    return float(q), nodeList
    

def getBestLouvainCommunities(graph, nIt = 1000, weighted=False):
    bestModularity = -float('inf')
    bestNodeList = None
    for i in range(nIt):
        if i%10 == 0:
            print "started ", i, "/", nIt
            sys.stdout.flush()
        q, nodeList = getLouvainCommunities(graph, weighted=weighted)
        if q > bestModularity:
            bestModularity=q
            bestNodeList = nodeList
    return {settings.modularity_abbr:bestModularity, settings.louvain_cluster_abbr:bestNodeList}

def computeClusterSimilarityMeasures(membershipLists, measures = settings.cluster_similarity_measures):
    """
    Computes pairwise clustering similarity measures between different clusterings (and conditions).
    
    Parameters:
        List of membership lists corresponding the clusterings
        
    Returns:
        A dict where the key corresponds to the similarity measure and value
        is a (upper triagonal) matrix containing the values of the similarity 
        measures.
    """
    membershipLists = np.array(membershipLists, dtype=int)
    nTot = len(membershipLists)
    resultDict = {}        
    for measure in measures:
        simMat = np.zeros( (nTot, nTot) )
        for i in range(0,nTot):
            clu1 = membershipLists[i]
            clu1 = [int(clu1[k]) for k in range(len(clu1))]
            for j in range(i,nTot):
                clu2 = membershipLists[j]
                clu2 = [int(clu2[k]) for k in range(len(clu2))]
                simMat[i,j] = igraph.compare_communities(clu1, clu2, measure)
                simMat[j,i] = simMat[i,j]
        resultDict[measure] = simMat
    return resultDict    
    
    