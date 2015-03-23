# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 15:52:03 2013

@author: rmkujala


A module for exporting data for other (mainly visualization uses).

Not in use / DEPRECATED!
"""
import settings
import dataio
import numpy as np
import scipy.io


# Examples for just one use case, not to be considered generic functions (yet)

def outPutSignificantLinksToGMLFlatNets():
    # load statistical data

    corrFDRStats = dataio.loadPickle(
        dataio.getCorrFDRStatsFileName())[settings.correlation_tag]

    corrFDRThresholdIndex = 8  # TODO: 'arbitrary, should be parametrizd
    print corrFDRStats[settings.fdrs_tag][corrFDRThresholdIndex], corrFDRStats[settings.thresholds_tag][corrFDRThresholdIndex]

    tValThresh = corrFDRStats[settings.thresholds_tag][corrFDRThresholdIndex]

    #dataio.loadPickle(dataio.getCorrPValStatsFileName() )
    if 0:  # change also to usage of qM
        qvalFDRStats = dataio.loadPickle(dataio.getCorrQValStatsFileName())
        qThreshold = 0.06  # TODO:  arbitrary, should be given as a parameter
        qMask = (qvalFDRStats[settings.qval_tag] < qThreshold)

    # not using qMask here
    treatmentMask = (corrFDRStats[settings.tval_tag] > tValThresh)  # *qMask
    controlMask = (corrFDRStats[settings.tval_tag] < -tValThresh)  # *qMask

    # load node loc data
    voxData = scipy.io.loadmat(
        params.fmriVoxelInfoFileName, squeeze_me=True)['rois']
    coords = voxData['flatmap']
    coords = np.array(list(coords))  # hack to make it one numpy array
    cosuce = voxData['cosuce']
    lr = voxData['lr']
    offset1nodes = (cosuce == 1) * (lr == 1)
    offset2nodes = (cosuce == 1) * (lr == 2)
    offset1 = np.array([-900, 50])
    offset2 = np.array([100, 50])
    # broadcasting
    coords += (offset1nodes.reshape(1, -1) * offset1.reshape(2, 1)).T
    coords += (offset2nodes.reshape(1, -1) * offset2.reshape(2, 1)).T

    okNodes = dataio.getOkIndices()
    coords = coords[okNodes]  # remove blacklist nodes
    cosuce = cosuce[okNodes]
    lr = lr[okNodes]
    nNodes = len(okNodes)

    linkIndexToNodeIndices = np.array(np.triu_indices(nNodes, 1)).T
    treatmentLinks = linkIndexToNodeIndices[treatmentMask.nonzero()[0]]
    controlLinks = linkIndexToNodeIndices[controlMask.nonzero()[0]]
    mode2Links = [treatmentLinks, controlLinks]

    # add extra marker nodes for visualization (after everything else)
    markers = np.array([[-900, 0], [900, 950]])
    # Using the background Enrico gave
    coords = np.vstack((coords, markers))
    # coords[:,1] *= -1    #at least when cytoscape is used

    # which links are outputted! TODO (better)!
    def acceptInterHemispehereLink(source, target):
        if cosuce[source] != 1 or cosuce[target] != 1:
            return False
        if lr[source] != lr[target]:
            return True
        return False

    def acceptIntraHemisphereLink(source, target):
        if cosuce[source] != 1 or cosuce[target] != 1:
            return False
        if lr[source] == lr[target]:
            return True
        return False

    def acceptMidHemisphereLink(source, target):
        if cosuce[source] == 1 and cosuce[target] == 1:
            return False
        return True

    for i, mode in enumerate(params.modes):
        rejectionLinks = mode2Links[i]

        # output links between hemispeheres + cortices
        interLinks = []
        intraLinks = []
        midLinks = []
        for link in rejectionLinks:
            source = link[0]
            target = link[1]
            if acceptInterHemispehereLink(source, target):
                interLinks.append(link)
            if acceptIntraHemisphereLink(source, target):
                intraLinks.append(link)
            if acceptMidHemisphereLink(source, target):
                midLinks.append(link)

        fname = params.outputDir + mode + "_inter.gml"
        outPutGML(interLinks, coords, fname)
        fname = params.outputDir + mode + "_intra.gml"
        outPutGML(intraLinks, coords, fname)
        fname = params.outputDir + mode + "_mid.gml"
        outPutGML(midLinks, coords, fname)


def outPutGML(links, nodeCoords2D, outFileName):
    """
    Outputs an adjacency matrix to gml format.
    Nodecoordinates are used for the layout.
    """
    f = open(outFileName, "w")
    strToWrite = 'graph [\nhierarchic 0\ndirected 0\nlabel "fMRI 2D-graph"\n'
    # write starters
    for i, nc in enumerate(nodeCoords2D):
        toAdd = "node [\n"
        toAdd += "id " + str(i) + "\n"
        toAdd += "graphics\t[\n"
        toAdd += "x " + str(nc[0]) + "\n"
        toAdd += "y " + str(nc[1]) + "\n"
        toAdd += "]\n"
        toAdd += "]\n"
        strToWrite += toAdd
        if i % 100 == 0:
            f.write(strToWrite)
            strToWrite = ""
    f.write(strToWrite)
    strToWrite = ""
    # write nodes
    # get edges as edge list
    # edges is a edge list
    for edge in links:
        toAdd = "edge [\n"
        toAdd += "source " + str(edge[0]) + "\n"
        toAdd += "target " + str(edge[1]) + "\n"
        toAdd += "]\n"
        strToWrite += toAdd
        if i % 100 == 0:
            f.write(strToWrite)
            strToWrite = ""
    strToWrite += "]"
    f.write(strToWrite)
    f.close()


def output_node_prop_to_nii(prop="strength", mask=True):
    propData = dataio.loadPickle(dataio.getNodeFDRStatsFileName())[prop]
    tVals = propData['tvalues']
    i = 3
    thresh = propData['thresh'][i]
    if mask:
        masking = ~(np.abs(tVals) < thresh)
        tVals = tVals * masking
    print thresh, propData['pfdrs'][i], propData['fdrs'][i], np.sum(np.abs(tVals) >= thresh)
    if mask:
        dataio.output_node_prop_to_nii(tVals, params.outputDir + prop + ".nii")
    else:
        dataio.output_node_prop_to_nii(
            tVals, params.outputDir + prop + "_unmasked.nii")


def outPutLouvainClustersToNii():
    for mode in params.modes:
        for i in range(1, 1 + params.modeToN[mode]):
            d = dataio.loadIndLouvainProperties(mode, i)
            data = d['louvain_clusters']
            data -= 0.5 * np.max(data)
            dataio.output_node_prop_to_nii(
                data, params.outputDir + "louvain_" + str(mode) + "_" + str(i) + ".nii")


def comstructuresToMapsWithLabelings(coms1, coms2, mapFileNames, fmriRoiInfoFileName, blacklistFileName):
    rois = dataio.loadMat(fmriRoiInfoFileName, squeeze_me=True)["rois"]
    okNodes = dataio.get_ok_nodes(len(rois), blacklistFileName)
    okRois = rois[okNodes]
    # weird empty array as first element...
    okRoiAalLabels = np.unique(okRois['aal_label'])[1:]
    okRoiAalIds = np.unique(okRois['aal_ID'])

    comIds1 = np.unique(coms1[okNodes])
    comIds2 = np.unique(coms2[okNodes])
    print comIds1
    overlapMatrix = np.zeros((len(comIds1), len(comIds2)))
    overlapNodesMatrix = []
    for i, comId1 in enumerate(comIds1):
        overlapNodesMatrix.append([])
        comId1Nodes = (coms1 == comId1)
        for j, comId2 in enumerate(comIds2):
            overlapNo
            comId2Nodes = (coms2 == comId2)
            # get joint nodes:
            jointNodes = comId1Nodes * comId2Nodes
            overlapNodesMatrix[i][j].append(np.nonzero(jointNodes))
            overlapMatrix[i, j] = np.sum(jointNodes)

    print np.sum(overlapMatrix)

#    print okRoiAalLabels, okRoiAalIds
    return None


def comstructureToMap(comStructure, mapFileName):
    """
    Takes in a (filtered) comstructure and outputs a .map file used for the generation fo the alluvial diagram.
    Args:
            comStructruce: nodewise cluster association
            mapFileName: the output file name
    """
    nNodes = len(comStructure)
    moduleLabels = np.sort(np.unique(comStructure))
    nModules = len(moduleLabels)
    modules = []
    moduleSizes = []
    for label in moduleLabels:
        modules.append(np.nonzero(comStructure == label)[0])
        moduleSizes.append(len(modules[-1]))
    outputString = "# modules: " + str(nModules) + "\n"
    outputString += "# modulelinks: 0\n"
    outputString += "# nodes: " + str(nNodes) + "\n"
    outputString += "# links: 0\n"
    outputString += "# codelength: 1.0\n"
    outputString += "*Directed \n"
    outputString += "*Modules " + str(nModules) + "\n"
    for i, moduleLabel in enumerate(moduleLabels):
        outputString += str(i + 1) + ' "Module' + str(i + 1) + \
            '"' + " " + str(moduleSizes[i]) + " 1.0\n"
    outputString += "*Nodes " + str(nNodes) + "\n"
    for i, module in enumerate(modules):
        for j, nodeLabelInt in enumerate(module):
            outputString += str(i + 1) + ":" + str(j + 1) + \
                ' "Node' + str(nodeLabelInt + 1) + '" 1.0\n'
    outputString += "*Links 0\n"
    f = open(mapFileName, "w")
    f.write(outputString)
    f.close()
