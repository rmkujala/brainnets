#!/usr/bin/env python
import numpy as np
from mlabwrap import mlab
#mlab.path(mlab.path(), 'ClusterPack-V1.0/')

from rpy2.robjects import r
from rpy2.robjects import IntVector
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

def matchClustersHungarianAlgo(cluster1,cluster2):
    '''
    Match the clusters using Hungarian Algorithm
    '''
    src=r.source('community_structure.R')
    hungarianmatch=r["hungarianmatch"]
    cluster=hungarianmatch(IntVector(cluster1),IntVector(cluster2))
    cluster=map(int,cluster)
    return [c-1 for c in cluster] 

def getClusterEnsembles(group=1,dataType='Best'):
    '''
    Determine the ensemble cluster for patients and groups
    dataType: 'Best', 'Consensus'
    '''
    if group == 1 :
        # Controls
        subjects=range(1,14)
    elif group== 2 :
        # Patients 
        subjects=range(14,27)

    if dataType=='Best':
        dirName='BestModules/'
        filePrefix='bestModule_subject_'
        commFilePrefix='bestModule_ClusterEnsembles_subject_'
    elif dataType=='Consensus':
        dirName='ConsensusModules/'
        filePrefix='consensusModule_subject_'
        commFilePrefix='consensusModule_ClusterEnsembles_subject_'

    communities={}
    for subjectNumber in subjects:
        commInFileName=dirName+filePrefix+str(subjectNumber)+'.txt'
        communities[subjectNumber]=np.loadtxt(commInFileName,dtype=int)

    commOutFileName=dirName+commFilePrefix+str(subjects[0])+'_'+str(subjects[-1])+'.txt'
    commOutFile=open(commOutFileName,'w')

    for k in communities.iterkeys():
        nNodes=communities[k].shape[1]
        break

    pRange=range(1,51)
    for idx, pLinks in enumerate(pRange):
        print pLinks
        comms=np.zeros((len(subjects),nNodes),dtype=int)
        for idx1, subjectNumber in enumerate(subjects):
            comms[idx1,:]=communities[subjectNumber][idx,:]
        
        bestPartition=mlab.mcla(comms+1)
        bestPartition=bestPartition-1
        nodeList=list(bestPartition.astype(int).flatten())
        nodeList=compressModuleNumber(nodeList)

        commOutFile.write(' '.join(map(str,nodeList))+'\n')
        commOutFile.flush()
    commOutFile.close()
    return commOutFileName

def getClusterEnsemblesBestModules():
    '''
    Determine the ensemble cluster for patients and groups
    Consider only the best partitions
    '''
    commOutFileName1=getClusterEnsembles(1,'Best')
    commOutFileName2=getClusterEnsembles(2,'Best')

    clusters1=np.loadtxt(commOutFileName1,dtype=int)
    clusters2=np.loadtxt(commOutFileName2,dtype=int)
    commOutFileName=commOutFileName2.replace('.txt','_Matched.txt')
    with open(commOutFileName,'w') as commOutFile:
        for i in range(clusters1.shape[0]):
            print i
            nodeList=matchClustersHungarianAlgo(clusters1[i,:]+1,clusters2[i,:]+1)
            commOutFile.write(' '.join(map(str,nodeList))+'\n')
            commOutFile.flush()

def getClusterEnsemblesConsensus():
    '''
    Determine the ensemble cluster for patients and groups
    Consider the consensus partitions
    '''
    commOutFileName1=getClusterEnsembles(1,'Consensus')
    commOutFileName2=getClusterEnsembles(2,'Consensus')

    clusters1=np.loadtxt(commOutFileName1,dtype=int)
    clusters2=np.loadtxt(commOutFileName2,dtype=int)
    commOutFileName=commOutFileName2.replace('.txt','_Matched.txt')
    with open(commOutFileName,'w') as commOutFile:
        for i in range(clusters1.shape[0]):
            print i
            nodeList=matchClustersHungarianAlgo(clusters1[i,:]+1,clusters2[i,:]+1)
            commOutFile.write(' '.join(map(str,nodeList))+'\n')
            commOutFile.flush()

def compressModuleNumber(nodeList):
    modifiedList=[]
    node_map = {}
    for node in nodeList:
        modifiedList.append(node_map.setdefault(node, len(node_map)))
    return modifiedList

def getGroupConsistency(dataType='Best'):
    '''
    Determine the consistence between the ensemble cluster for patients and groups
    dataType: 'Best', 'Consensus'
    '''
    if dataType=='Best':
        dirName='BestModules/'
        filePrefix='bestModule_subject_'
        commFilePrefix='bestModule_ClusterEnsembles_subject_'
        commOutFilePrefix='bestModule_ClusterEnsembles_consistency_subject_'
    elif dataType=='Consensus':
        dirName='ConsensusModules/'
        filePrefix='consensusModule_subject_'
        commFilePrefix='consensusModule_ClusterEnsembles_subject_'
        commOutFilePrefix='consensusModule_ClusterEnsembles_consistency_subject_'

    pRange=range(1,51)
    for group in [1,2]:
        if group == 1 :
            # Controls
            subjects=range(1,14)
        elif group== 2 :
            # Patients 
            subjects=range(14,27)

        communities={}
        for subjectNumber in subjects:
            commInFileName=dirName+filePrefix+str(subjectNumber)+'.txt'
            communities[subjectNumber]=np.loadtxt(commInFileName,dtype=int)

        commInFileName=dirName+commFilePrefix+str(subjects[0])+'_'+str(subjects[-1])+'.txt'
        commonCommunities=np.loadtxt(commInFileName,dtype=int)

        commOutFileName=dirName+commOutFilePrefix+str(subjects[0])+'_'+str(subjects[-1])+'.txt'

        nNodes=commonCommunities.shape[1]

        with open(commOutFileName,'w') as commOutFile:
            for idx, pLinks in enumerate(pRange):
                print pLinks
                comms=np.zeros(nNodes,dtype=int)
                for idx1, subjectNumber in enumerate(subjects):
                    nodeList=matchClustersHungarianAlgo(commonCommunities[idx,:]+1,communities[subjectNumber][idx,:]+1)
                    nodeList=np.array(nodeList)
                    comms += (commonCommunities[idx,:]==nodeList).astype(int)
                commOutFile.write(' '.join(map(str,comms))+'\n')
                commOutFile.flush()


    subjects=range(1,14)
    commInFileName=dirName+commFilePrefix+str(subjects[0])+'_'+str(subjects[-1])+'.txt'
    comm1=np.loadtxt(commInFileName,dtype=int)
    subjects=range(14,27)
    commInFileName=dirName+commFilePrefix+str(subjects[0])+'_'+str(subjects[-1])+'.txt'
    commInFileName=commInFileName.replace('.txt','_Matched.txt')
    comm2=np.loadtxt(commInFileName,dtype=int)
    commOutFileName=dirName+commOutFilePrefix.replace('consistency_subject_','difference')+'.txt'
    np.savetxt(commOutFileName,(comm1!=comm2).astype(int),fmt='%d')

def getGroupConsistencyAndDifferenceBestModules():
    '''
    Determine the consistence and differences between the ensemble cluster for patients and groups
    Consider the best partitions
    '''
    getGroupConsistency('Best')

def getGroupConsistencyAndDifferenceConsensusModules():
    '''
    Determine the consistence and differences between the ensemble cluster for patients and groups
    Consider the consensus partitions
    '''
    getGroupConsistency('Consensus')

def getGroupConsistencySignificance(group=1,dataType='Best',nIterations=1000):
    '''
    Determine the nodes that are significatly consistent in the ensemble cluster for patients or controls
    dataType: 'Best', 'Consensus'
    '''
    if dataType=='Best':
        dirName='BestModules/'
        filePrefix='bestModule_subject_'
        commFilePrefix='bestModule_ClusterEnsembles_subject_'
        commOutFilePrefix='bestModule_ClusterEnsembles_consistencySignificance_subject_'
    elif dataType=='Consensus':
        dirName='ConsensusModules/'
        filePrefix='consensusModule_subject_'
        commFilePrefix='consensusModule_ClusterEnsembles_subject_'
        commOutFilePrefix='consensusModule_ClusterEnsembles_consistencySignificance_subject_'

    if group == 1 :
        # Controls
        subjects=range(1,14)
    elif group== 2 :
        # Patients 
        subjects=range(14,27)

    pRange=range(35,51)
    communities={}
    for subjectNumber in subjects:
        commInFileName=dirName+filePrefix+str(subjectNumber)+'.txt'
        communities[subjectNumber]=np.loadtxt(commInFileName,dtype=int)

    commInFileName=dirName+commFilePrefix+str(subjects[0])+'_'+str(subjects[-1])+'.txt'
    commonCommunities=np.loadtxt(commInFileName,dtype=int)
    nNodes=commonCommunities.shape[1]

    commOutFileName=dirName+commOutFilePrefix+str(subjects[0])+'_'+str(subjects[-1])+'.txt'

    with open(commOutFileName,'a') as commOutFile:
        for idx, pLinks in enumerate(pRange):
            print pLinks
            significanceComms=np.zeros(nNodes)
            for nIter in range(nIterations):
                print nIter
                comms=np.zeros((len(subjects),nNodes),dtype=int)
                randSubjects=np.random.randint(subjects[0],subjects[-1]+1,size=len(subjects))
                for idx1, subjectNumber in enumerate(randSubjects):
                    comms[idx1,:]=communities[subjectNumber][idx,:]

                bestPartition=mlab.mcla(comms+1)
                bestPartition=bestPartition-1
                nodeList=list(bestPartition.astype(int).flatten())
                nodeList=np.array(compressModuleNumber(nodeList))

                nodeList=matchClustersHungarianAlgo(commonCommunities[idx,:]+1,nodeList+1)
                nodeList=np.array(nodeList)
                significanceComms += (commonCommunities[idx,:]==nodeList).astype(int)
            significanceComms /= nIterations
            commOutFile.write(' '.join("%0.4f" % x for x in significanceComms)+'\n')
            commOutFile.flush()

def getGroupDifferenceSignificance(dataType='Best',nIterations=1000):
    '''
    Determine if the difference between cluster ensamble are significant or not
    '''
    if dataType=='Best':
        dirName='BestModules/'
        filePrefix='bestModule_subject_'
        commOutFilePrefix='bestModule_ClusterEnsembles_differenceSignificance'
    elif dataType=='Consensus':
        dirName='ConsensusModules/'
        filePrefix='consensusModule_subject_'
        commOutFilePrefix='consensusModule_ClusterEnsembles_differenceSignificance'

    pRange=range(43,51)
    subjects=range(1,27)
    communities={}
    for subjectNumber in subjects:
        commInFileName=dirName+filePrefix+str(subjectNumber)+'.txt'
        communities[subjectNumber]=np.loadtxt(commInFileName,dtype=int)

    for k in communities.iterkeys():
        nNodes=communities[k].shape[1]
        break

    commOutFileName=dirName+commOutFilePrefix+'.txt'
    with open(commOutFileName,'a') as commOutFile:
        for idx, pLinks in enumerate(pRange):
            print pLinks
            significanceComms=np.zeros(nNodes)
            for nIter in range(nIterations):
                print nIter
                subjects=np.random.permutation(subjects)
                comms=np.zeros((len(subjects),nNodes),dtype=int)
                for idx1, subjectNumber in enumerate(subjects[:13]):
                    comms[idx1,:]=communities[subjectNumber][idx,:]

                bestPartition=mlab.mcla(comms+1)
                bestPartition=bestPartition-1
                nodeList=list(bestPartition.astype(int).flatten())
                nodeList1=np.array(compressModuleNumber(nodeList))

                comms=np.zeros((len(subjects),nNodes),dtype=int)
                for idx1, subjectNumber in enumerate(subjects[13:]):
                    comms[idx1,:]=communities[subjectNumber][idx,:]

                bestPartition=mlab.mcla(comms+1)
                bestPartition=bestPartition-1
                nodeList=list(bestPartition.astype(int).flatten())
                nodeList2=np.array(compressModuleNumber(nodeList))

                nodeList2=matchClustersHungarianAlgo(nodeList1+1,nodeList2+1)
                nodeList2=np.array(nodeList2)
                significanceComms += (nodeList1!=nodeList2).astype(int)
            significanceComms /= nIterations
            commOutFile.write(' '.join("%0.4f" % x for x in significanceComms)+'\n')
            commOutFile.flush()

def getGroupConsistencyAndDifferenceSignificanceBestModules():
    '''
    Determine the consistence and differences between the ensemble cluster for patients and groups
    Consider the best partitions
    '''
    getGroupConsistencySignificance(group=1,dataType='Best',nIterations=1000)
    getGroupConsistencySignificance(group=2,dataType='Best',nIterations=1000)
    getGroupDifferenceSignificance(dataType='Best',nIterations=1000)

def getGroupConsistencyAndDifferenceSignificanceConsensusModules():
    '''
    Determine the consistence and differences between the ensemble cluster for patients and groups
    Consider the best partitions
    '''
    getGroupConsistencySignificance(group=1,dataType='Consensus',nIterations=1000)
    getGroupConsistencySignificance(group=2,dataType='Consensus',nIterations=1000)
    getGroupDifferenceSignificance(dataType='Consensus',nIterations=1000)
