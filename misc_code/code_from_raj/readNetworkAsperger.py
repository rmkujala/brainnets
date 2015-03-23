import numpy as np
import scipy.io
import cPickle as pickle
from netpython import pynet, transforms

def nodeMapping():
    '''
    Get the node mapping
    '''
    fileName='/triton/becs/scratch/networks/rajkp/DataSets/AS_networks/nodeMapping.pck'
    try :
        with open(fileName) as f:
            nodeMap=pickle.load(f)
    except IOError :
        nNodes=6235
        blackList=scipy.io.loadmat('/scratch/braindata/jsalmi/AS-fMRI/corrected_nifti/data_motion_regressed/AS_blacklist.mat', squeeze_me=True, struct_as_record=False)
        blackList=blackList['blacklist']
        allNodes=np.ones(nNodes,dtype=bool)
        allNodes[blackList-1]=False
        
        nodeMap={}
        with open(fileName.replace('.pck','.txt'),'w') as fp:
            for i,node in enumerate((allNodes==True).nonzero()[0]):
                nodeMap[node+1]=i
                fp.write(('%d %d\n') % (node, i) )
        with open(fileName, 'w') as f:
                pickle.dump(nodeMap, f)
    return nodeMap


def nodeLeftRightHemisphereLabel():
    '''
    For each node check dtermine the whether it is in left or right hemisphere
    '''
    fileName='/triton/becs/scratch/networks/rajkp/DataSets/AS_networks/nodeLeftRightHemisphereLabel.pck'
    try :
        with open(fileName) as f:
            nodeLabel=pickle.load(f)
    except IOError :
        rois=scipy.io.loadmat('/triton/becs/archive/braindata/2011/networks/networks/rois_6mm_plus.mat', squeeze_me=True)
        labels=rois['rois']['label']
        nodeLabel={}
        for i,l in enumerate(labels):
            nodeLabel[i+1]=str(l[-1])
        with open(fileName, 'w') as f:
            pickle.dump(nodeLabel, f)
    return nodeLabel

def getEuclideanDistanceMatrix(nodeMap):
    '''
    Get the Euclidean Distance between the pair of nodes
    '''
    distArray=np.loadtxt('/scratch/networks/rajkp/DataSets/funpsy/euclideanDistance_roi_6mm.txt')
    nodeList=np.array(sorted(nodeMap.keys()))-1
    distArray=distArray[nodeList,:][:,nodeList]
    return distArray

def getNetwork(fileName,nodeMap):
    '''
    Read the network and return the link weights as an array
    '''
    inDirName='/triton/becs/scratch/networks/rajkp/DataSets/AS_networks/'
    net=pynet.SymmNet()
    with open(inDirName+fileName) as fp:
        for line in fp:
            n1,n2,w=line.split('\t')
            n1,n2,w=int(n1),int(n2),float(w)
            try :
                net[nodeMap[n1],nodeMap[n2]]=w
            except KeyError :
                pass
    mst=transforms.mst(net,maximum=True)
    netArray=np.zeros((len(net.edges)-len(mst.edges),), dtype=[('node1', 'i4'), ('node2', 'i4'), ('weight','f8')])

    counter=0
    for n1, n2, w in net.edges: 
        if not mst[n1,n2]:
            netArray[counter]=n1,n2,w
            counter += 1

    netArray=np.sort(netArray,order=['weight'])[::-1]
    counter=0
    mstArray=np.zeros((len(mst.edges),), dtype=[('node1', 'i4'), ('node2', 'i4'), ('weight','f8')])
    for n1, n2, w in mst.edges: 
        mstArray[counter]=n1,n2,w
        counter += 1

    return net, np.hstack((mstArray,netArray))

def getNetworkMatrix(fileName,nodeMap):
    '''
    Get the Euclidean Distance between the pair of nodes
    '''
    nNodes=6235
    inDirName='/triton/becs/scratch/braindata/hihalme/AS_networks/motion_regressed/original_nets/'
    mat = scipy.io.loadmat(inDirName+fileName)
    net=mat['adj']

    nodeList=np.array(sorted(nodeMap.keys()))-1
    net=net[nodeList,:][:,nodeList]
    return net
