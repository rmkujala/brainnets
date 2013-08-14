# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 10:08:05 2013

@author: rmkujala
"""
import dataio
import gencomps
import settings
import params
import numpy as np
#from scipy.stats import gaussian_kde
import statistics
import matplotlib.pyplot as plt
#import matplotlib as mpl
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import host_subplot
#import mpl_toolkits.axisartist as AA


from matplotlib import rc
rc('text', usetex=True)
#rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Sans serif']})
rc('font',**{'family':'serif','serif':['Computer Modern serif']})
rc('legend',fontsize=12)
#rc[]

def rerunAllPlots():
    """
    Helper function, which runs all the visualizations.
    """
    #use this if there has been any changes in colors etc.
    plotCorrelationDists()
    globalWPropsPlot()
    globalUWPropsPlot()
    plotLouvainSimilarityMeasures()

    
def getLinBins(nBins, binMin, binMax):
    bins = np.linspace(binMin, binMax, nBins+1)
    binCenters = 0.5*(bins[:-1]+bins[1:])
    return bins, binCenters

def getBinCounts(values, bins):
    return np.bincount(np.digitize(values, bins), minlength=len(bins))[1:]

def plotCorrelationDists():
    """
    Plot the (joint/combined) correlation distributions for movie and rest
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    nBins = 100
    corrBins, corrBinCenters= getLinBins(nBins,-1,1.)
    binCounts = {"movie":np.zeros(nBins), "rest":np.zeros(nBins)}

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, mode in enumerate(params.modes):
        for j in range(params.modeToN[mode]):
            print j
            flatCorrMat = dataio.loadFlattenedCorrMatrix(mode, j+1)
            binCounts[mode] += getBinCounts(flatCorrMat, corrBins)
        #normalize
        binCounts[mode] = binCounts[mode]/(np.sum(binCounts[mode]*(corrBins[1]-corrBins[0])))
        ax.plot(corrBinCenters, binCounts[mode], color = params.colors[i], label=mode)

    ax.set_xlabel(settings.getPropTexName("corr"))
    ax.set_ylabel(r"Probability density P(c)")
    ax.legend(loc=0)
    plt.savefig(params.outputDir+"corrDists.pdf", format="pdf")
    
def plotCorrTValDist():
    """
    Plot the (joint/combined) correlation distributions for movie and rest
    """
    #get tvals
    corrStats = dataio.loadPickle(dataio.getCorrFDRStatsFilename(paired=True, asym=False))
    tVals = corrStats[settings.tval_abbr]
    minVal = np.min(tVals)
    maxVal = np.max(tVals)
    nBins = 100
    bins, binCenters= getLinBins(nBins,minVal-(np.abs(minVal)*0.1), maxVal+(np.abs(maxVal)*0.1))
    binCounts = getBinCounts(tVals, bins)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #normalize
    binCounts = binCounts*1./(np.sum(binCounts)*(bins[1]-bins[0]))
    ax.plot(binCenters, binCounts, label=r""+params.mode1+"-"+params.mode2)
    ax.set_xlabel(settings.getPropTexName(settings.tval_abbr)
    ax.set_ylabel(r"Probability density")
    ax.legend(loc=0)
    plt.savefig(params.outputDir+"tvalDist.pdf", format="pdf")
    
#def plotInvCDF(ax, values, xscale='lin', yscale='log', nSteps=1000):
#    """
#    UNDER CONSTRUCTION
#    """
#    print "plotInvCDF under construction, returning...."
#    return
#    stepPoints = None
#    if nBins == 'all':
#        #warning: you do not want to plot like 10**8 points!
#        nBins = len(values)
#        stepPoints = np.sort(values)
#    
#    values.sort()
#    values = nBins

def plotLinProfile(ax, values, nBins=100, color="r"):
    """
    Plots the 
    PDF and CDF in linear space (x)
    """
    values = np.array(values)
    valMax = np.max(values)
    valMin = np.min(values)
    linBins = np.linspace(valMin-(0.1*np.abs(valMin)), valMax+(0.1*np.abs(valMax)), nBins+1)
    bins, binCenters = getLinBins(nBins, valMin-0.1*np.abs(valMin), valMax+0.1*np.abs(valMax) ) 
    binCounts = np.bincount(np.digitize(values, linBins), minlength=nBins+1  )[1:]
    #normalize
    binCounts = binCounts*1./(len(values)*(linBins[1]-linBins[0]))
    ax.plot(binCenters, binCounts, color=color)


def getGlobalPropDataToPlot(weighted=False):
    """
    Returns:
            propMeans: a python dictionary with keys as the property
            -> propMeans[prop][mode][p] = mean
            propMeanConfIntervals
            -> propData[prop][mode][p] = [low, high]
            propStats: a python dict containing the perm. ttest pvalues
            -> propStats[prop][p] = pvalue
            percentages: a list of the percentages
    """
    if weighted:
        data = dataio.mergeAndLoadGlobalWProperties()
        statsData = dataio.loadPickle(dataio.getGlobalWStatsFilename())
    else:
        data = dataio.mergeAndLoadGlobalUWProperties()
        statsData = dataio.loadPickle(dataio.getGlobalUWStatsFilename())
    
    propMeans = {}
    propMeanConfIntervals = {}
    propPVals = {}
    percentages = data[settings.percentages_abbr]
    
    for prop in data.keys():
        propData = data[prop]
        propStats = statsData[prop]
        if prop != settings.percentages_abbr:
            modeMeanDict = {}
            modeMeanDict[params.modes[0]] = np.average(propData[:params.n1,:], 0)
            modeMeanDict[params.modes[1]] = np.average(propData[params.n1:,:], 0)
            propMeans[prop] = modeMeanDict

            modeMeanErrDict = {}
            modeMeanErrDict[params.modes[0]] = statistics.bootsrapMeanConfInterval( propData[:params.n1,:].T, percent=params.bootStrapCoverageInPercents)
            modeMeanErrDict[params.modes[1]] = statistics.bootsrapMeanConfInterval( propData[params.n1:,:].T, percent=params.bootStrapCoverageInPercents)            
            propMeanConfIntervals[prop] = modeMeanErrDict
            propPVals[prop] = propStats[settings.pval_abbr]
    return propMeans, propMeanConfIntervals, propPVals, percentages

def globalPropsPlot(data, filenamePostFix, x_len=None, y_len = None):
    propMeans, propMeanConfIntervals, propPVals, percentages = data
    nProps = len(propMeans.keys())
    if x_len == None:
        x_len = int(np.sqrt(nProps)+1)
    if y_len == None:
        y_len = int(np.ceil(nProps/x_len))
    fig = plt.figure(figsize = (x_len*4,y_len*4) )
    color0 = params.colors[0]
    color1 = params.colors[1]
    for i, key in enumerate(propMeans):
        ax = host_subplot(y_len, x_len, i)
        means0 = propMeans[key][params.mode1]
        means1 = propMeans[key][params.mode2]
        confInt0 = propMeanConfIntervals[key][params.mode1].T
        confInt1 = propMeanConfIntervals[key][params.mode2].T
        _plotComparisonAndPValue(ax, percentages, means0, means1, confInt0, confInt1, color0, color1, propPVals[key], key )
    plt.tight_layout(pad=2.2)
    fig.savefig(params.outputDir+filenamePostFix, format="pdf")
    
def globalWPropsPlot():
    data = getGlobalPropDataToPlot(weighted=True)
    globalPropsPlot(data, "globalWProps.pdf",x_len=4, y_len=1)
    
def globalUWPropsPlot():
    data = getGlobalPropDataToPlot(weighted=False)
    globalPropsPlot(data, "globalUWProps.pdf")
    
    
def _plotComparisonAndPValue(ax, percentages, means0, means1, confInt0, confInt1, color0, color1, pVals, measure_key, xlabels = []):
    pax = ax.twinx()
    
    ax.set_xlabel(settings.getPropTexName(settings.percentages_abbr))
    ax.set_ylabel(settings.getPropTexName(measure_key))
    #pax.set_ylabel("P-value")
    
    p1 = ax.plot(percentages, means0, label=settings.mrModes[0], color = color0, lw=1.5)
    p2 = ax.plot(percentages, means1, label=settings.mrModes[1], color = color1, lw=1.5)
    alpha = 0.25
    ax.fill_between(percentages, confInt0[0], confInt0[1], color=color0, alpha=alpha)
    ax.fill_between(percentages, confInt1[0], confInt1[1], color=color1, alpha=alpha) 
    p3 = pax.semilogy(percentages, pVals, "-o", color="0.35", markersize = 2.5, label="p-value", zorder=-10)
    if params.studyType == "paired":
        pax.axhline(y=2**-(params.n1-1), xmin=0, xmax=1, color='0.6', lw=1.5, ls="--")
    pax.set_ylim( (10**-4, 10**0) )
    pax.yaxis.tick_right() #to take the left side log scale off (a bug in matplotlib with loglog)
    
    # make the legend
    lns = p1+p2+p3
    labs = [l.get_label() for l in lns]
    l = ax.legend(lns, labs, loc=0, numpoints=1,handlelength = 1)
    l.get_frame().set_alpha(0.75)
    
def plotLouvainSimilarityMeasures():
    data = dataio.mergeAndLoadLouvainProperties()
    clusterings = data[settings.louvain_cluster_abbr]
    cluster_n = []
    for i, c in enumerate(clusterings):
        cluster_n.append(len(np.unique(c)))
            
    print cluster_n
    #print statistics.permutationTest(statistics.meanDifference, data[settings.modularity_abbr], 13, paired='all')
    #print statistics.permutationTest(statistics.meanDifference, np.array(cluster_n), 13, paired='all')
        
    simMatricesDict = gencomps.computeClusterSimilarityMeasures(clusterings)
        
    for i,measure in enumerate(simMatricesDict):
        fig = plt.figure(figsize=(8,10))
        ax = fig.add_subplot(1, 1,1)
        plotSimMatrix(ax, simMatricesDict[measure], measure)
        commonX = 0.1
        commonY = 0.93
        controlX = 0.69
        controlY = 0.48        
        treatmentX = 0.34
        treatmentY = 0.76
        fig.text(treatmentX,commonY, params.mode1TexName, rotation='horizontal', va='center',  ha='center')
        fig.text(commonX,treatmentY, params.mode1TexName, rotation='vertical', va='center',  ha='center')
        fig.text(controlX,commonY, params.mode2TexName, rotation='horizontal', va='center',  ha='center')
        fig.text(commonX,controlY, params.mode2TexName, rotation='vertical', va='center',  ha='center')
        
        print "permuation test statistics for the similarity measure " + measure
        stats = statistics.groupMeanDiffSimMatrixPermutationTest(simMatricesDict[measure], params.n1, True)
        
        means, meanConfIntervals = statistics.bootstrapMeanConfIntervalForSimMatrix(simMatricesDict[measure], settings.mrNsubjs)
        treatmentText = r"\begin!center?mean = {:.3f} \\ 95\%: ({:.3f}-{:.3f})\end!center?".format(means[0], meanConfIntervals[0][0], meanConfIntervals[0][1])
        controlText = r"\begin!center?mean = {:.3f} \\ 95\%: ({:.3f}-{:.3f})\end!center?".format(means[1], meanConfIntervals[1][0], meanConfIntervals[1][1])
        treatmentText = treatmentText.replace("!", "{")
        controlText = controlText.replace("!", "{")
        treatmentText = treatmentText.replace("?", "}")
        controlText = controlText.replace("?", "}")
        
        fig.text(treatmentX-0.05, treatmentY-0.05, treatmentText, bbox={"facecolor":"w", "alpha":0.8}, va='center',  ha='center', alpha=0.8)
        fig.text(controlX-0.05, controlY-0.05, controlText, bbox={"facecolor":"w", "alpha":0.5}, va='center',  ha='center')
        fig.text(0.50, 0.30, r"permutation test: p = {:.5f}".format(stats[1][2]), va='center',  ha='center')
        fig.savefig(params.outputDir+"louvain_"+measure+".pdf", format="pdf")

    
def plotSimMatrix(ax, m, measure):
    """
    Plot a similarity matrix in to the given axis "ax". "m" contains the image
    and "measure" is the name of the measure to be plotted (ie. a string)
    """
    sortedvals = np.sort(np.unique(m.flatten()))
    if measure == 'vi':
        cmap = cm.hot
        vmax = np.max(m)
        vmin = sortedvals[1]
    elif measure == 'nmi' or measure == 'adjusted_rand':
        cmap = cm.hot_r
        sortedvals = np.sort(np.unique(m.flatten()))
        vmax_start = sortedvals[-1]
        for i in range(2,len(sortedvals)):
            vmax = sortedvals[-i]
            if (vmax_start-vmax)/vmax_start > 0.0001:
                break
        vmin = np.min(m)
    else:
        print "trying to plot unknown measure with function _plotSimMatrix..."
        
    im = ax.imshow(m, interpolation='nearest', cmap=cmap, vmax =vmax, vmin=vmin)
    #cbarax = mpl.colorbar.make_axes(ax)
    #mpl.colorbar.Colorbar(ax, im)
    #, cax=cbarax)
    #cbax = mpl.colorbar.make_axes(parent, **kw)¶
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal')
    
    cbar.set_label(settings.getPropTexName(measure))
    ticks = np.arange(0, len(m))
    ax.yaxis.set_ticklabels(ticks+1)
    ax.yaxis.set_ticks(ticks)
    ax.xaxis.set_ticks(ticks)
    ax.xaxis.set_ticklabels(ticks+1)
    
    ax.xaxis.set_ticks_position('top')
    ax.tick_params(length=0, width=0, colors='k')
    ymed = np.average(ax.get_ylim())
    xmed = np.average(ax.get_xlim())
    ax.axhline(y=ymed, xmin=0, xmax=1, color='0.0', lw=1.0, ls="-")
    ax.axvline(x=xmed, ymin=0, ymax=1, color="0.0", lw=1.0, ls="-")
