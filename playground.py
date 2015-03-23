"""
This module contains code, which is not currently in real use or is just being tested.
"""


class TreeNode(object):

    """ For plotting dendrograms, not fully functional and unused"""

    def __init__(self, simTresh, children=[], parent=None):
        self.simTresh = simTresh
        self.children = children
        self.parent = parent
        self.highest_known_parent = parent

    def addParent(self, parent):
        self.parent = parent
        self.highest_known_parent = parent

    def getRoot(self):
        if self.highest_known_parent == None:
            return self
        else:
            self.highest_known_parent = self.highest_known_parent.getRoot()
            return self.highest_known_parent

    def getLeafCount(self):
        if self.children:
            leafs = 0
            for child in self.children:
                leafs += child.getLeafCount()
            return leafs
        else:
            return 1

    def plotDendrogram(self, ax, xRange, lowestSim, maxSim):
        if not self.children:
            xmid = np.average(xRange)
            ax.plot([xmid, xmid], [self.simTresh, lowestSim])
            return xmid
        else:
            yCur = self.simTresh
            child1, child2 = self.children
            child1Leafs = child1.getLeafCount()
            child2Leafs = child2.getLeafCount()
            child1xfrac = float(child1Leafs) / (child1Leafs + child2Leafs)
#            print child1Leafs, child2Leafs
            child1xRange = [
                xRange[0], xRange[0] + (xRange[1] - xRange[0]) * child1xfrac]
            child2xRange = [
                xRange[0] + (xRange[1] - xRange[0]) * child1xfrac, xRange[1]]
            #child1x = np.average(child1xRange)
            xmid1 = child1.plotDendrogram(ax, child1xRange, lowestSim, maxSim)
            xmid2 = child2.plotDendrogram(ax, child2xRange, lowestSim, maxSim)
            ax.plot([xmid1, xmid2], [yCur, yCur], "k")
            ax.plot([xmid1, xmid1], [yCur, child1.simTresh], "k")
            ax.plot([xmid2, xmid2], [yCur, child2.simTresh], "k")
            return np.average([xmid1, xmid2])


def makeClusteringTree(adjMat):
    """ For plotting dendrograms, not fully functional """
    import sys
    nNodes = len(adjMat)
    sys.setrecursionlimit(np.maximum(nNodes + 10, 10000))
    allNodes = []
    for i in range(nNodes):
        allNodes.append(TreeNode(1.0))
    triu_indices = np.triu_indices_from(adjMat, 1)
    triu_vals = adjMat[triu_indices]
    sortedtriuval_args = np.argsort(triu_vals)[::-1]
    triu_indices = np.array(triu_indices)
    nJoins = 0
    for arg in sortedtriuval_args:
        i, j = triu_indices[:, arg]
        nodeI = allNodes[i]
        nodeJ = allNodes[j]
        nodeiRoot = nodeI.getRoot()
        nodejRoot = nodeJ.getRoot()
        if nodeiRoot == nodejRoot:
            continue
        else:
            val = triu_vals[arg]
            newRoot = TreeNode(val, children=[nodeiRoot, nodejRoot])
            nodeiRoot.addParent(newRoot)
            nodejRoot.addParent(newRoot)
            nJoins += 1
            if nJoins == nNodes - 1:
                return newRoot


def plotDendrogramFromRoot(rootNode, ax=None):
    """ For plotting dendrograms, not fully functional """
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    xRange = np.array([0, 1])
    rootNode.plotDendrogram(ax, xRange, 1, 0)
    return ax.get_figure(), ax


def plotDendrogram(adjMat):
    """ For plotting dendrograms, not fully functional """
    root = makeClusteringTree(adjMat)
    return plotDendrogramFromRoot(root)
