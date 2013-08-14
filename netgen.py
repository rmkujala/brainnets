import igraph
import numpy as np

def makeNetWithLargestLinksAndMST(corrMat, percentage=2.0, weighted=False):
    """
    Creates a network taking certain percentage of links having the 
    largest weights based on the given correlation/adjacency matrix.
    If weighted = False creates an unweighted graph
    If weighted = True creates an weighted graph with w_{ij} = 1/\rho_{ij}
    """
    m = corrMat.copy()
    n = len(m)
    initEdges = []    
    weights = []
    for i in range(0,n): 
        for j in range(i+1,n):
            initEdges.append((i,j))
            weights.append(m[i,j])            
    weights = np.array(weights)
    g = igraph.Graph(n, initEdges, directed=False)
    mst = g.spanning_tree(-1*weights, return_tree=False)
    #pick up the minimum spanning tree (-1* due that above)   
    edgelist = []
    for ei in mst:
        edge = g.es[ei]
        source = edge.source
        target = edge.target        
        edgelist.append( (source, target) )
        m[source, target] = -1.0 #take these links away from the orig. mat
    nLinksMax = (n*(n-1))/2    
    #how many links we still need to take:
    nLinksYetToTake = np.max([int(nLinksMax*percentage/100.0-(n-1)),0]) # mst already there
    mflat = m.flatten()
    #flatten for finding largest elements
    flatindices = mflat.argsort()[::-1][:nLinksYetToTake]
    for flatindex in flatindices:
        edgelist.append(np.unravel_index(flatindex, (n,n)))        

    if weighted == False:    
        #make the new _binary_ graph
        binGraph = igraph.Graph(n)
        binGraph.add_edges(edgelist)
        return binGraph    
    if weighted == True:
        graph = igraph.Graph(n)
        graph.add_edges(edgelist)
        weights = [corrMat[edge[0],edge[1]] for edge in edgelist]
        graph.es["weight"] = weights
        return graph    
