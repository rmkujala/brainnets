#
# This module contains the study specific parameters and settings
# which should be edited before running any scripts.
#
#
# EDITING THESE VARIABLES INSIDE A SCRIPT CAN BE HAZARDOUS!
# (at least for the time being)

import os

#STUDY SETUP
##########################################################################################

#number of subjects in the 1st group (treatment)
n1 = None #e.g. 13 (an int)
mode1 = None #change to suitable, e.g. "movie" or "treatment" 
#the names of the correlation matrices should have the name
# inputDir+mode1+"_"+str(i)+".mat"

#number of subjects in the 2nd group (control, if existing)
n2 = None #e.g. 13 (an int)
mode2 = None #change to suitable, e.g. "rest" or "control" 
modes = [mode1, mode2]
#It can get hazardous if n1 or n2 is changed in a python script!:
modeToN = {mode1:n1, mode2:n2}
#densities of the graphs to be built in percentages
pRange = range(1,31)
pSpecial = 2.

#STATISTICAL SETUP (used for permutation testing):
# are the subjects matched or not?
pairedStudy = True/False
if pairedStudy == True:
    assert n1==n2, "With paired setup you should have equal number of subjects in both groups"

permSamples = 1e5 #option permSamples = 'all' is also availble if 'paired' setting is used
if permSamples == 'all':
    assert pairedStudy, "if option 'all' is specified for the nIt, you should have a 'paired' studyType"
    
#RUN PARAMETERS
louvain_iterations = 1000
parallel_cpus = 12 #good for triton

#DATAINFO
##########################################################################################
#number of voxels initially
nVoxels = -1 #number of 
fmriVoxelInfoFileName = None #e.g. "/proj/braindata/eglerean/net_viz/rois6mm_final_plus_IF2.mat"
blacklistNodes = None # e.g. "blacklist.mat" 

#PLOTTING AND VISUALISATION:
##########################################################################################
color1 = 'r'
color2 = 'b'
colors = [color1, color2]
mode1TexName = r""+str(mode1) #change to something more appropriate, if wanted
mode2TexName = r""+str(mode2) #change to something more appropriate, if wanted
bootStrapCoverageInPercents = 95 
bootStrapSamples = 1e6 #(usually 1e6 still not too expensive)

#DIRECTORIES:
##########################################################################################
#name of the directory where everything should be outputted
outputDir = os.getcwd()+"/"
#name of the directory where all data should be loaded from
inputDir = os.getcwd()+"/"
#name of the directory for temporary computation files
tmpDir = "/tmp/"

