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

#################
#Set directories#
#################
today = str(date.today())
#read1 = raw_input('Read 1 fastq Location (/dir/R1.fastq): ')
#read2 = raw_input('Read 2 fastq Location (/dir/R2.fastq): ')
#outputDir = raw_input('Output Location (/dir): ')
infoOutput = raw_input('Info Writeup About This Run (info/n): ')

#just hardcoding to expedite testing
read1 = '/home/alex/Dropbox/Code/FERMI/testInput/R1.fastq'
read2 = '/home/alex/Dropbox/Code/FERMI/testInput/R2.fastq'
outputDir = '/home/alex/Dropbox/Code/FERMI/testOutput'
#read1 = '/media/alex/Extra/Dropbox/Code/FERMI/testInput/R1.fastq'
#read2 = '/media/alex/Extra/Dropbox/Code/FERMI/testInput/R2.fastq'
#outputDir = '/media/alex/Extra/Dropbox/Code/FERMI/testOutput'
previousDict = raw_input('Would you like to load previously sorted data? (Y/n): ')
if previousDict == 'Y':
    prevDictLoc = raw_input('Location of previous sorted data (/dir/data.pkl): ')

if not path.exists(outputDir):
    #make the output directory expanduser is used to allow ~/Desktop shortcuts
    mkdir(os.path.expanduser(outputDir))

####################
#Set Run Parameters#
####################
useDefaults = raw_input('Use Default Parameters? (Y/n): ')

if useDefaults == 'Y':
    #number of mismatches allowed when calling two UMIs the same
    distance_stringency = 1
    #threshold % of reads that must contain a given base read
    varThresh = 0.75
    #the number of required supporting reads of each UMI pair
    supportingReads = 5

elif useDefaults == 'n':
    distance_stringency = int(raw_input('Allowed UMI Mismatches (1): '))
    varThresh = float(raw_input('Read Prevalence Threshold (0.75): '))
    supportingReads = int(raw_input('Required Supporting Reads (5): '))

#create output directory
outputDir = outputDir + '/' + today + '_' + str(supportingReads) + '_' + str(varThresh)

clusterRun = raw_input('Submit Run to LRS Cluster? (Y/n): ')

########################
#Write Dated Output Dir#
########################
#this now allows for multiple daily runs with the same parameters
if path.exists(outputDir):
    counter = 1
    tempDir = outputDir
    while path.exists(tempDir):
        tempDir = outputDir
        tempDir = tempDir + '_Run_' + str(counter)
        counter += 1
    outputDir = outputDir + '_Run_' + str(counter - 1)

if not path.exists(outputDir):
    #make the output directory
    mkdir(os.path.expanduser(outputDir))

twoUmiOut = outputDir + '/twoUMIs.fastq'
final_output_file = outputDir + '/finalOutput.fastq'
coverage_file = outputDir + '/coverageData.txt'
infoFile = outputDir + '/runInfo.txt'
parametersUsed = outputDir + '/parametersUsed.txt'
pickleOutput = outputDir + '/sortedSeqData.pkl'

#############################
#Record Files and Parameters#
#############################
if infoOutput != 'n':
    info = open(infoFile, 'w')
    info.write(infoOutput)
    info.close()

target = open(parametersUsed, 'w')
target.write("Read 1 Location: %s\n" %(read1))
target.write("Read 2 Location: %s\n" %(read2))
target.write("Distance Stringency: %d\n" %(distance_stringency))
target.write("Variant Threshold: %f\n" %(varThresh))
target.write("Supporting Reads: %d\n" %(supportingReads))
target.close()

##############################
#Output All Variables to File#
##############################

#################
#Run Main Script#
#################
if clusterRun == 'Y':
    pass
elif clusterRun == 'n':
    #without the following check main.py gets run forever
    if __name__ == "__main__":
        system("python main.py")
