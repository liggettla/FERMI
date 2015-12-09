from datetime import date
import os
from os import mkdir
from os import path
from concatenateUMI import concatenateUMI
from goodCollapseDictionary import buildNestedDict
from goodCollapseDictionary import collapseNestedDict
import pickle
from outputCoverage import outputCov
import pdb
from os import system


from goodCollapseDictionary import buildListDict
from goodCollapseDictionary import collapseReadsListDict

#####################
#Concatate R1 and R2#
#####################
#attach 3' UMI from R2 onto R1 read
#this is a necessary step to process reads with 100-150 cycle chemistry
#200 cycle chemistry makes this unnecessary
concatenateUMI(read1, read2, twoUmiOut)

##############################
#Build/Get Seq Data Structure#
##############################
#build dict binning reads by concatenated UMIs
if previousDict == 'n':
    seqDict = buildNestedDict(twoUmiOut, distance_stringency, pickleOutput)
#retrieve previous seq data structure
elif previousDict == 'Y':
    prevData = open(prevDictLoc, 'rb')
    seqDict = pickle.load(prevData)
    prevData.close()

################
#Collapse Reads#
################
#calculate individual read length between the UMIs
with open(twoUmiOut, 'r') as target:
    header = next(target)
    readSeq = next(target).rstrip('\n')
    readLength = len(readSeq) - 12

#collapse reads on binned UMI data structure
collapseNestedDict(seqDict, varThresh, final_output_file, supportingReads, readLength)

#####################
#Output Seq Coverage#
#####################
outputCov(twoUmiOut, final_output_file, distance_stringency, coverage_file)

###################
#Plot Allele Freqs#
###################
#this should output a plot automatically, but is not runtime tested
#does not output log10() of plot yet
if False:
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    r = ro.r
    r.setwd(outputDir)
    f = r('read.table("allelefreqs.txt", header = FALSE)')
    grdevices = importr('grDevices')
    grdevices.png(file="alleleFreq.png", width=800, height=500)
    r.hist(f[0], breaks=100, main = '5 Reads', xlab='Variant Freq', ylab='# Vars')
    grdevices.dev_off()

###################
#Run Shell Scripts#
###################
'''
This is not yet function, but shows the general idea of how to run shell scripts
autonomously and still pass them necessary variables.
The idea is that variables are output to some variable.txt file like so:
    output
    /dir/here
    read1
    /dir/R1.fastq
    read2
    /dir/R2.fastq
Then the following code should actually be run from within the actual shell scripts
and it will pull out necessary variables through a combined grep/tail search
'''
if False:
    read1=$(cat variables.txt | grep -A 1 'read1' variables.txt | tail -n 1)
