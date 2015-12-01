from datetime import date
import os
from os import mkdir
from concatenateUMI import concatenateUMI
from goodCollapseDictionary import buildListDict
from goodCollapseDictionary import collapseReadsListDict
from outputCoverage import outputCov
import pdb

###############
#Home Computer#
###############
#read1 = '/home/alex/Desktop/testData/A2_small1.fastq'
#read2 = '/home/alex/Desktop/testData/A2_small2.fastq'
#twoUmiOut = '/home/alex/Desktop/testData/twoUMIs.fastq'
#final_output_file = '/home/alex/Desktop/testData/finalOutput.fastq'
#coverage_file = '/home/alex/Desktop/testData/coverageData'

##############
#Lab Computer#
##############
#read1 = '/home/alex/Desktop/testData/A2_small1.fastq' #R1 from paired end reads
#read2 = '/home/alex/Desktop/testData/A2_small2.fastq' #R2 from paired end reads
#twoUmiOut = '/home/alex/Desktop/testData/twoUMIs.fastq' #File generated by combining R1 and R2
#final_output_file = '/home/alex/Desktop/testData/finalOutput.fastq'
#coverage_file = '/home/alex/Desktop/testData/coverageData'

#########
#Cluster#
#########
#read1 = '/vol3/home/liggettl/TruSeqPanel/8.17.2015_HiSeqFastq/unzipped/A2_GAAACCC_L008_R1_001.fastq'
#read2 = '/vol3/home/liggettl/TruSeqPanel/8.17.2015_HiSeqFastq/unzipped/A2_GAAACCC_L008_R2_001.fastq'
#twoUmiOut = '/vol3/home/liggettl/TruSeqPanel/8.17.2015_HiSeqFastq/5_Reads_75_Percent/twoUMIs.fastq'
#final_output_file = '/vol3/home/liggettl/TruSeqPanel/8.17.2015_HiSeqFastq/10_Reads_99_Percent/finalOutput.fastq'
#coverage_file = '/vol3/home/liggettl/TruSeqPanel/8.17.2015_HiSeqFastq/10_Reads_99_Percent/coverageData'

#################
#Set directories#
#################

today = str(date.today())
read1 = raw_input('Read 1 fastq Location (/dir/R1.fastq): ')
read2 = raw_input('Read 2 fastq Location (/dir/R2.fastq): ')
outputDir = raw_input('Output Location (/dir): ')
infoOutput = raw_input('Info Writeup About This Run (info/n): ')

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
    distance_stringency = raw_input('Allowed UMI Mismatches (1): ')
    varThresh = raw_input('Read Prevalence Threshold (0.75): ')
    supportingReads = raw_input('Required Supporting Reads (5): ')

#create output directory
outputDir = outputDir + '/' + today + '_' + str(supportingReads) + '_' + str(varThresh)
#make the output directory
mkdir(os.path.expanduser(outputDir))

twoUmiOut = outputDir + '/twoUMIs.fastq'
final_output_file = outputDir + '/finalOutput.fastq'
coverage_file = outputDir + '/coverageData.txt'
infoFile = outputDir + '/runInfo.txt'
parametersUsed = outputDir + '/parametersUsed.txt'

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

#attach 3' UMI from R2 onto R1 read
#this is a necessary step to process reads with 100-150 cycle chemistry
#200 cycle chemistry makes this unnecessary
concatenateUMI(read1, read2, twoUmiOut)

#build dict binning reads by concatenated UMIs
seqDict = buildListDict(twoUmiOut, distance_stringency)

#calculate individual read length between the UMIs
with open(twoUmiOut, 'r') as target:
    header = next(target)
    readSeq = next(target).rstrip('\n')
    readLength = len(readSeq) - 12

#collapse reads on binned UMIs
collapseReadsListDict(seqDict, varThresh, final_output_file, supportingReads, readLength)

#calculate capture coverage
outputCov(twoUmiOut, final_output_file, distance_stringency, coverage_file)
