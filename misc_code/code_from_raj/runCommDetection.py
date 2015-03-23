#!/usr/bin/env python
import numpy as np
import scipy.io
import cPickle as pickle
from netpython import pynet, transforms
from mlabwrap import mlab
#mlab.path(mlab.path(), 'ClusterPack-V1.0/')

import os
import time
import subprocess as sub
from generalFn import frange
import time
from netpython import netio
import modulesProperties
import itertools
from collections import defaultdict
from readNetworkAsperger import nodeMapping, getEuclideanDistanceMatrix, getNetwork

def generateThresholdedNet(netArray,nLinks,subjectNumber=1):
    '''
    Generate the thresholded network and the binary file 
    netArray = Network as an array
    nLinks = Number of Links
    '''
    fileName='graph_'+str(subjectNumber)+'.edg'
    tNet=pynet.SymmNet()
    with open(fileName,'w') as fp:
        for n1, n2, w in netArray[:nLinks]:
            fp.write(("%d %d %f\n") % (n1,n2,w))
            tNet[n1,n2]=w

    runstr='./convert -i '+fileName+' -o '+fileName.replace('.edg','.bin')+' -w '+fileName.replace('.edg','.weights')
    p = sub.call(runstr.split())
    return tNet

def getAverageLinkWeight(netArray,nLinks):
    '''
    Get the average link weight
    '''
    return netArray['weight'][:nLinks].mean()

def detectCommunities(commFileName,weighted=False):
    '''
    Determine the modules at the final hierarchical level
    '''
    # Detect the communities
    if weighted==True :
        runstr = './community '+commFileName.replace('.tree','.bin')+' -l -1 -r 1 -w '+commFileName.replace('.tree','.weights')
    else :
        runstr = './community '+commFileName.replace('.tree','.bin')+' -l -1 -r 1'
    with open(commFileName,'w') as fp:
        p = sub.Popen(runstr.split(),stdout=fp,stderr=sub.PIPE)
    _,q=p.communicate()

    # Get the number of hierarchical level
    runstr = './hierarchy ' + commFileName
    p = sub.Popen(runstr.split(),stdout=sub.PIPE,stderr=sub.PIPE)
    hLevels,_=p.communicate()
    hLevel=int(hLevels.split('\n')[0].split(':')[1])-1
    # Nodes in the final hierarchical level
    runstr = './hierarchy ' + commFileName + ' -l ' + str(hLevel) 
    p = sub.Popen(runstr.split(),stdout=sub.PIPE,stderr=sub.PIPE)
    nodeOutput,_=p.communicate()
    
    # Read the community 
    nodeList=[]
    for line in nodeOutput.split('\n')[:-1]: 
        nodeList.append(int(line.split()[1]))

    return float(q), nodeList

def getModularitySpecificFraction(subjectNumber,nIterations=1000):
    '''
    More detailed modularity comparision for some specific values
    '''
    dirName='ConsensusModules/'
    outFileName=dirName+'avgModuleProps_subject_'+str(subjectNumber)+'.txt'
    outFile=open(outFileName,'a')

    print 'Subject: ', subjectNumber
    pRange=[1,2,5,10,15,20]

    nodeMap=nodeMapping()
    fileName=str(subjectNumber)+'.edg'
    net,netArray=getNetwork(fileName,nodeMap)
    nNodes=len(net)
    distArray=getEuclideanDistanceMatrix(nodeMap)
    outFile.write('#fracOfLinks numberOfLinks avgLinkWeight modularity+-std numberOfModules+-std avgCompactness+-std nullModelProbability+-std\n')

    for idx, pLinks in enumerate(pRange):
        outFileName1=dirName+'allModules_subject_'+str(subjectNumber)+'_pLinks_'+str(pLinks)+'.txt'
        outFile1=open(outFileName1,'w')
        print '\npLinks: ', pLinks
        nLinks=(pLinks*nNodes*(nNodes-1))/(2*100)
        outFile.write(('%d %d ') % (pLinks,nLinks))
        # Generate the binarry network
        tNet=generateThresholdedNet(netArray,nLinks,subjectNumber)
        avgLinkWt=getAverageLinkWeight(netArray,nLinks)
        outFile.write(('%0.4f ') % (avgLinkWt))
        consensusMatrix=np.zeros((nNodes,nNodes),dtype=int)

        mProps=np.zeros((nIterations,), dtype=[('modularity', 'f4'), ('nModules', 'i4'), ('nInLinks','i4'),('dInLinks','f4'),('avgWtInLinks','f4'),('avgCompactness','f4'),('threshold','f4')])
        commFileName='graph_'+fileName.replace('.edg','.tree')
        print 'nIterations:',
        for nIter in range(nIterations):
            if nIter%10 == 0:
                print '%d' % nIter
            q, nodeList = detectCommunities(commFileName)

            outFile1.write(' '.join(map(str,nodeList))+'\n')
            outFile1.flush()

            nComms=max(nodeList)+1
            commProps=np.zeros((nComms,), dtype=[('size','i4')])
            commProps['size']=np.bincount(nodeList)

            mProps[nIter]['modularity'] = q 

            # Number of communities
            mProps[nIter]['nModules']=commProps.shape[0]

            # Get the consensus matrix
            consensusMatrix=modulesProperties.getConsensusMatrix(nodeList,consensusMatrix)
            # Consensus matrix expectation (Null model)
            mProps[nIter]['threshold']=modulesProperties.getThreholdConsensusMatrix(commProps)

            #return nodeList, distArray
            commCompactness=modulesProperties.getModularCompactness(nodeList,distArray)
            mProps[nIter]['avgCompactness']=modulesProperties.getAverageModularCompactness(commCompactness,commProps)

        for prop in ['modularity', 'nModules', 'avgCompactness','threshold']:
            outFile.write(('%0.4f %0.4f ') % (mProps[prop].mean(),mProps[prop].std()))
        outFile.write('\n')
        outFile.flush()

        commOutFileName=dirName+'consensusMatrix_subject_'+str(subjectNumber)+'_pLinks_'+str(pLinks)+'.txt'
        np.savetxt(commOutFileName,consensusMatrix,fmt='%d')
        outFile1.close()
    outFile.close()

def getModularityComparision(subjectNumber,nIterations=100):
    '''
    Determine the best modules
    Compare different properties of these modules
    Calculate :
    - modularity, 
    - number of modules, 
    - number of inter-modular link,
    - weighted density of inter-modular links, 
    - Avg weight of internal links
    - Avg compactness
    '''

    dirName='BestModules/'
    outFileName=dirName+'bestModuleProps_subject_'+str(subjectNumber)+'.txt'
    commOutFileName=dirName+'bestModule_subject_'+str(subjectNumber)+'.txt'
    outFile=open(outFileName,'w')
    commOutFile=open(commOutFileName,'w')

    print 'Subject: ', subjectNumber
    pRange=range(1,51)
    nodeMap=nodeMapping()
    fileName=str(subjectNumber)+'.edg'
    net,netArray=getNetwork(fileName,nodeMap)
    nNodes=len(net)
    distArray=getEuclideanDistanceMatrix(nodeMap)
    outFile.write('#fracOfLinks numberOfLinks avgLinkWeight modularity numberOfModules numberOfInternalLinks densityOfInternalLinks avgWeightInternalLinks avgCompactness\n')
    for idx, pLinks in enumerate(pRange):
        print 'pLinks: ', pLinks
        nLinks=(pLinks*nNodes*(nNodes-1))/(2*100)
        outFile.write(('%d %d ') % (pLinks,nLinks))
        # Generate the binarry network
        tNet=generateThresholdedNet(netArray,nLinks,subjectNumber)
        avgLinkWt=getAverageLinkWeight(netArray,nLinks)
        outFile.write(('%0.4f ') % (avgLinkWt))
        commFileName='graph_'+fileName.replace('.edg','.tree')
        bestQ=None
        print 'nIterations:',
        for nIter in range(nIterations):
            q, nodeList = detectCommunities(commFileName)

            #return tNet, nodeList
            if bestQ==None or bestQ<q :
                bestQ=q
                bestPartition=nodeList
            if nIter%10 == 0:
                print '(%d %0.4f)' % (nIter,q)

        print ''
        # Calculate only for best module
        nodeList=bestPartition
        # Get various properties of the modules
        commProps=modulesProperties.communityProperties(tNet,nodeList)
        # Modularity of the partition
        modularity = modulesProperties.getModularity(commProps, nLinks)
        # Number of internal links 
        nInLinks = modulesProperties.getTotalInLinks(commProps)
        # Density of internal links
        dInLinks = modulesProperties.getDensityInLinks(commProps)
        # Average weight of internal links
        avgWtInLinks=modulesProperties.getInLinksAvgWeight(commProps)
        # Number of communities
        nModules=commProps.shape[0]
        # Average Compactness
        commCompactness=modulesProperties.getModularCompactness(nodeList,distArray)
        avgCompactness=modulesProperties.getAverageModularCompactness(commCompactness,commProps)
        outFile.write(('%0.4f %d %d %0.4f %0.4f %0.4f\n') % (modularity,nModules,nInLinks,dInLinks,avgWtInLinks,avgCompactness))
        commOutFile.write(' '.join(map(str,nodeList))+'\n')
        outFile.flush()
        commOutFile.flush()
    outFile.close()
    commOutFile.close()

def getModularityComparisionConsensus(subjectNumber,nIterations=100):
    '''
    Determine the best modules
    Compare different properties of these modules
    Calculate :
    - modularity, 
    - number of modules, 
    - number of inter-modular link,
    - weighted density of inter-modular links, 
    - Avg weight of internal links
    - Avg compactness
    '''

    dirName='ConsensusModules/'
    outFileName=dirName+'consensusModuleProps_subject_'+str(subjectNumber)+'.txt'
    commOutFileName=dirName+'consensusModule_subject_'+str(subjectNumber)+'.txt'
    outFile=open(outFileName,'w')
    commOutFile=open(commOutFileName,'w')

    print 'Subject: ', subjectNumber
    pRange=range(1,51)
    nodeMap=nodeMapping()
    fileName=str(subjectNumber)+'.edg'
    net,netArray=getNetwork(fileName,nodeMap)
    nNodes=len(net)
    distArray=getEuclideanDistanceMatrix(nodeMap)
    outFile.write('#fracOfLinks numberOfLinks avgLinkWeight modularity numberOfModules numberOfInternalLinks densityOfInternalLinks avgWeightInternalLinks avgCompactness\n')
    for idx, pLinks in enumerate(pRange):
        print 'pLinks: ', pLinks
        nLinks=(pLinks*nNodes*(nNodes-1))/(2*100)
        outFile.write(('%d %d ') % (pLinks,nLinks))
        # Generate the binarry network
        tNet=generateThresholdedNet(netArray,nLinks,subjectNumber)
        avgLinkWt=getAverageLinkWeight(netArray,nLinks)
        outFile.write(('%0.4f ') % (avgLinkWt))
        commFileName='graph_'+fileName.replace('.edg','.tree')
        bestQ=None
        nodeListArray=np.zeros((nIterations,nNodes),dtype=int)
        print 'nIterations:',
        for nIter in range(nIterations):
            q, nodeList = detectCommunities(commFileName)

            nodeListArray[nIter,:]=nodeList
            if bestQ==None or bestQ<q :
                bestQ=q
                bestPartition=nodeList
            if nIter%10 == 0:
                print '(%d %0.4f)' % (nIter,q)

        print ''
        bestPartition=mlab.mcla(nodeListArray+1)
        bestPartition=bestPartition-1
        # Calculate only for best module
        nodeList=list(bestPartition.astype(int).flatten())

        nodeList=compressModuleNumber(nodeList)
        # Get various properties of the modules
        commProps=modulesProperties.communityProperties(tNet,nodeList)
        # Modularity of the partition
        modularity = modulesProperties.getModularity(commProps, nLinks)
        # Number of internal links 
        nInLinks = modulesProperties.getTotalInLinks(commProps)
        # Density of internal links
        dInLinks = modulesProperties.getDensityInLinks(commProps)
        # Average weight of internal links
        avgWtInLinks=modulesProperties.getInLinksAvgWeight(commProps)
        # Number of communities
        nModules=commProps.shape[0]
        # Average Compactness
        commCompactness=modulesProperties.getModularCompactness(nodeList,distArray)
        avgCompactness=modulesProperties.getAverageModularCompactness(commCompactness,commProps)
        outFile.write(('%0.4f %d %d %0.4f %0.4f %0.4f\n') % (modularity,nModules,nInLinks,dInLinks,avgWtInLinks,avgCompactness))
        commOutFile.write(' '.join(map(str,nodeList))+'\n')
        outFile.flush()
        commOutFile.flush()
    outFile.close()
    commOutFile.close()

def compressModuleNumber(nodeList):
    modifiedList=[]
    node_map = {}
    for node in nodeList:
        modifiedList.append(node_map.setdefault(node, len(node_map)))
    return modifiedList

#def getConsensusModules(net,netArray,subjectNumber,nIterations=100):
def getConsensusModules(subjectNumber,nIterations=100):
    '''
    Get the consensus modules
    Check if the modules are block diagonal
    '''
    nodeMap=nodeMapping()
    fileName=str(subjectNumber)+'.edg'
    print "Reading the network"
    net,netArray=getNetwork(fileName,nodeMap)
    nNodes=len(net)
    print "Creating Dist matrix"
    distArray=getEuclideanDistanceMatrix(nodeMap)

    dirName='ConsensusModules/'
    inFileName=dirName+'avgModuleProps_subject_'+str(subjectNumber)+'.txt'
    commOutFileName=dirName+'consensusModule_subject_'+str(subjectNumber)+'.txt'
    outFileName=dirName+'consensusModuleProps_subject_'+str(subjectNumber)+'.txt'
    outFile=open(outFileName,'w')
    commOutFile=open(commOutFileName,'w')

    print "Creating consensus network"
    outFile.write('#fracOfLinks numberOfLinks avgLinkWeight modularity numberOfModules numberOfInternalLinks densityOfInternalLinks avgWeightInternalLinks avgCompactness\n')
    randProb={}
    with open(inFileName) as fp:
        fp.readline()
        for line in fp:
            line1=line.split()
            randProb[int(line1[0])]=float(line1[9])

    pRange=[1,2,5,10,15,20]
    for pLinks in pRange:
        print '\npLinks: ', pLinks
        nLinks=(pLinks*nNodes*(nNodes-1))/(2*100)
        outFile.write(('%d %d ') % (pLinks,nLinks))
        tNet=generateThresholdedNet(netArray,nLinks,subjectNumber)
        avgLinkWt=getAverageLinkWeight(netArray,nLinks)
        outFile.write(('%0.4f ') % (avgLinkWt))

        commInFileName=dirName+'consensusMatrix_subject_'+str(subjectNumber)+'_pLinks_'+str(pLinks)+'.txt'
        getNetworkFromTxtFile(commInFileName,randProb[pLinks],subjectNumber)
        commFileName='graph_'+fileName.replace('.edg','.tree')
        bestQ=None
        for nIter in range(nIterations):
            q, nodeList=detectCommunities(commFileName,weighted=True)
            if bestQ==None or bestQ<q :
                bestQ=q
                bestPartition=nodeList
            if nIter%10 == 0:
                print '(%d %0.4f)' % (nIter,q)

        print ''
        # Calculate only for best module
        nodeList=bestPartition

        # Get various properties of the modules
        commProps=modulesProperties.communityProperties(tNet,nodeList)
        # Modularity of the partition
        modularity = modulesProperties.getModularity(commProps, nLinks)
        # Number of internal links 
        nInLinks = modulesProperties.getTotalInLinks(commProps)
        # Density of internal links
        dInLinks = modulesProperties.getDensityInLinks(commProps)
        # Average weight of internal links
        avgWtInLinks=modulesProperties.getInLinksAvgWeight(commProps)
        # Number of communities
        nModules=commProps.shape[0]
        # Average Compactness
        commCompactness=modulesProperties.getModularCompactness(nodeList,distArray)
        avgCompactness=modulesProperties.getAverageModularCompactness(commCompactness,commProps)
        outFile.write(('%0.4f %d %d %0.4f %0.4f %0.4f\n') % (modularity,nModules,nInLinks,dInLinks,avgWtInLinks,avgCompactness))
        commOutFile.write(' '.join(map(str,nodeList))+'\n')
        outFile.flush()
        commOutFile.flush()
    outFile.close()
    commOutFile.close()

def getNetworkFromTxtFile(fileName,threshold=0,subjectNumber=1,loops=False):
    '''
    Generate the network from the txt-file (adjacency list)
    threshold : ignore link weights lower than threshold
    Remove the nan and assuming that the matrix is symmetric
    loops = False (ignore self loops); True (include self loops)
    Returns an undirected graph
    '''
    mat=np.loadtxt(fileName)
    mat=mat/1000.
    idx=(mat>=threshold).nonzero()
    fileName='graph_'+str(subjectNumber)+'.edg'
    with open(fileName,'w') as fp:
        for i, j in itertools.izip(idx[0],idx[1]):
            fp.write(("%d %d %f\n") % (i,j,mat[i,j]))

    runstr='./convert -i '+fileName+' -o '+fileName.replace('.edg','.bin')+' -w '+fileName.replace('.edg','.weights')
    p = sub.call(runstr.split())

def getNetworkFromTxtFileGroup(pLinks,subjects,randProb,loops=False):
    '''
    Generate the network from the txt-file (adjacency list)
    threshold : ignore link weights lower than threshold
    Remove the nan and assuming that the matrix is symmetric
    loops = False (ignore self loops); True (include self loops)
    Returns an undirected graph
    For the group
    '''
    thresholds=[randProb[s][pLinks] for s in subjects]
    threshold=np.mean(thresholds) 
    matGroup=None
    dirName='ConsensusModules/'
    for subjectNumber in subjects:
        print subjectNumber
        fileName=dirName+'consensusMatrix_subject_'+str(subjectNumber)+'_pLinks_'+str(pLinks)+'.txt'
        mat=np.loadtxt(fileName)
        if matGroup==None:
            matGroup=mat
        else:
            matGroup+= mat

    mat=matGroup/(len(subjects)*1000.)
    idx=(matGroup>=threshold).nonzero()
    fileName='graph_'+str(subjects[0])+'_'+str(subjects[-1])+'.edg'
    with open(fileName,'w') as fp:
        for i, j in itertools.izip(idx[0],idx[1]):
            fp.write(("%d %d %f\n") % (i,j,mat[i,j]))

    runstr='./convert -i '+fileName+' -o '+fileName.replace('.edg','.bin')+' -w '+fileName.replace('.edg','.weights')
    p = sub.call(runstr.split())

def getConsensusModulesGroups(group,nIterations=100):
    '''
    Get the consensus modules
    Check if the modules are block diagonal
    '''
    if group == 1 :
        subjects=range(1,14)
    elif group== 2 :
        subjects=range(14,27)

    nodeMap=nodeMapping()
    print "Creating Dist matrix"
    distArray=getEuclideanDistanceMatrix(nodeMap)

    dirName='ConsensusModules/'
    randProb=defaultdict(dict)
    for subjectNumber in subjects:
        inFileName=dirName+'avgModuleProps_subject_'+str(subjectNumber)+'.txt'
        with open(inFileName) as fp:
            fp.readline()
            for line in fp:
                line1=line.split()
                randProb[subjectNumber][int(line1[0])]=float(line1[9])

    commOutFileName=dirName+'consensusModule_subject_'+str(subjects[0])+'_'+str(subjects[-1])+'.txt'
    outFileName=dirName+'consensusModuleProps_subject_'+str(subjects[0])+'_'+str(subjects[-1])+'.txt'
    outFile=open(outFileName,'w')
    commOutFile=open(commOutFileName,'w')

    print "Creating consensus network"
    outFile.write('#fracOfLinks numberOfLinks avgLinkWeight modularity numberOfModules numberOfInternalLinks densityOfInternalLinks avgWeightInternalLinks avgCompactness\n')

    pRange=[1,2,5,10,15,20]
    nodeListAll={}
    for pLinks in pRange:
        print '\npLinks: ', pLinks

        getNetworkFromTxtFileGroup(pLinks,subjects,randProb)
        commFileName='graph_'+str(subjects[0])+'_'+str(subjects[-1])+'.tree'
        bestQ=None
        for nIter in range(nIterations):
            q, nodeList=detectCommunities(commFileName,weighted=True)
            if bestQ==None or bestQ<q :
                bestQ=q
                bestPartition=nodeList
            if nIter%10 == 0:
                print '(%d %0.4f)' % (nIter,q)

        print ''
        # Calculate only for best module
        nodeListAll[pLinks]=bestPartition
        commOutFile.write(' '.join(map(str,nodeListAll[pLinks]))+'\n')
        commOutFile.flush()
    commOutFile.close()

    avgLinkWt,modularity,nModules,nInLinks,dInLinks,avgWtInLinks,avgCompactness=defaultdict(list),defaultdict(list),defaultdict(list),defaultdict(list),defaultdict(list),defaultdict(list),defaultdict(list)
    for subjectNumber in subjects:
        fileName=str(subjectNumber)+'.edg'
        print "Reading the network"
        net,netArray=getNetwork(fileName,nodeMap)
        nNodes=len(net)
        for pLinks in pRange:
            nLinks=(pLinks*nNodes*(nNodes-1))/(2*100)
            tNet=generateThresholdedNet(netArray,nLinks,subjectNumber)
            avgLinkWt[pLinks].append(getAverageLinkWeight(netArray,nLinks))

            # Get various properties of the modules
            commProps=modulesProperties.communityProperties(tNet,nodeListAll[pLinks])
            # Modularity of the partition
            modularity[pLinks].append(modulesProperties.getModularity(commProps, nLinks))
            # Number of internal links 
            nInLinks[pLinks].append(modulesProperties.getTotalInLinks(commProps))
            # Density of internal links
            dInLinks[pLinks].append(modulesProperties.getDensityInLinks(commProps))
            # Average weight of internal links
            avgWtInLinks[pLinks].append(modulesProperties.getInLinksAvgWeight(commProps))
            # Number of communities
            nModules[pLinks].append(commProps.shape[0])
            # Average Compactness
            commCompactness=modulesProperties.getModularCompactness(nodeList,distArray)
            avgCompactness[pLinks].append(modulesProperties.getAverageModularCompactness(commCompactness,commProps))

    for pLinks in pRange:
        nLinks=(pLinks*nNodes*(nNodes-1))/(2*100)
        outFile.write(('%d %d ') % (pLinks,nLinks))
        avgL,m,nM,nInL,dInL,avgWtInL,avgC=[np.mean(x[pLinks]) for x in [avgLinkWt,modularity,nModules,nInLinks,dInLinks,avgWtInLinks,avgCompactness]]
        outFile.write(('%0.4f %0.4f %d %d %0.4f %0.4f %0.4f\n') % (avgL,m,nM,nInL,dInL,avgWtInL,avgC))
        outFile.flush()
    outFile.close()
