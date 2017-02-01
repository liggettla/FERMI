#!/usr/bin/env python
'''
Ensure that both of the following have been run on cluster headnode if
clusterRun is going to be utilized:
module load bwa
module load vt
'''

from datetime import date
import os
from os import mkdir
from os import path
#from concatenateUMI import concatenateUMI
#from goodCollapseDictionary import buildNestedDict
#from goodCollapseDictionary import collapseNestedDict
import pickle
#from outputCoverage import outputCov
#import pdb
from os import system
import os.path
import time

#from goodCollapseDictionary import buildListDict
#from goodCollapseDictionary import collapseReadsListDict

#####################
#Argparse Flag Input#
#####################
import argparse

#parser.add_argument('--info', '-i', required=True, type=int, help='specifies the readlimit')
#parser.add_argument('--something', '-s', action='store_true', help='no')

parser = argparse.ArgumentParser()
parser.add_argument('--nfo', '-n', type=str, help='Info writeup about a particular run that will be output in the run directory.')
parser.add_argument('--largefiles', '-l', action='store_true', help='Outputs all generated fastq files generated during analysis.')
parser.add_argument('--avoidalign', '-a', action='store_true', help='Only runs through initial analysis of input fastq files, and does not align to reference or call variants.')
parser.add_argument('--outdir', '-o', required=True, type=str, help='Specifies output directory where all analysis files will be written')
parser.add_argument('--indir', '-i', required=True, type=str, help='Specifies the input directory that contains the fastq files to be analyzed.')
parser.add_argument('--single', '-s', action='store_true', help='Only process a single set of paired end reads.')
parser.add_argument('--prevdict', '-p', type=str, help='Specify a previously output pickle file containing collapsed fastq data as an input instead of raw fastq files.')
parser.add_argument('--umimismatch', '-u', type=int, help='Specify the number of mismatches allowed in a UMI pair to still consider as the same UMI')
parser.add_argument('--varthresh', '-v', type=float, help='Specify the percentage of reads that must contain a particular base for that base to be used in the final consensus read')
parser.add_argument('--readsupport', '-r', type=int, help='Specifies the number of reads that must have a given UMI sequence in order to be binned as a true capture event, and not be thrown out.')
parser.add_argument('--clustersubmit', '-c', action='store_true', help='Submit run to cluster computing rather than running locally')
parser.add_argument('--filterao', '-f', type=int, help='Specifies the AO cuttoff for reported variants, where -f 5 would eliminate all variants that are seen 5 times or less. Default == 5.')
parser.add_argument('--dpfilter', '-d', type=int, help='Read depth elimination threshold. If specified as -d 500 only variants found in a locus read greater than 500 times will be reported. Default == 500.')
parser.add_argument('--freebayes', '-b', type=str, help='Location of freebayes in the format of /dir/freebayes')
parser.add_argument('--errorrate', '-e', action='store_true', help='Overall pcr amplification + sequencing error rates will be estimated and returned')
parser.add_argument('--readLength', '-q', type=int, help='Manually set the read length. If this is not set, length will be automatically set as the number of bases found between the two UMI sequences.')
parser.add_argument('--badBaseSubstitute', '-x', action='store_true', help='This flag will trigger replacement of bad bases with N instead of invalidating an entire capture.')
parser.add_argument('--reference', '-y', type=str, help='Set the location of the human reference genome hg19.fa and supporting files.')

args = parser.parse_args()

#################
#Set directories#
#################
# info writeup
today = str(date.today())
if args.nfo:
    infoOutput = args.nfo
else:
    infoOutput = 'n'
#infoOutput = raw_input('Info Writeup About This Run (info/n): ')
#infoOutput = 'n'

# large fastq file writing
if args.largefiles:
    noBigFiles = 'n'
else:
    noBigFiles = 'Y'
#noBigFiles = raw_input('Suppress Output of Large Files (Y/n): ')

# avoid alignment and var calling
if args.avoidalign:
    alignAndVar = 'n'
else:
    alignAndVar = 'Y'
#alignAndVar = raw_input('Align and Call Variants (Y/n): ')

# set output directory
outputDir = args.outdir
#outputDir = raw_input('Output Location (/dir): ')
#outputDir = './testOutput'

# reference genome location
REF = args.reference

###################
#Get Input File(s)#
##################
if args.single:
    numFiles = 'Y'
else:
    numFiles = 'n'
#numFiles = raw_input('Process only single file? (Y/n): ')

# set input directory
inputDir = args.indir
#inputDir = raw_input('Input Dir (/dir): ')

readList={}

if numFiles == 'Y':
    read1 = raw_input('Read 1 fastq Location (R1.fastq): ')
    read2 = raw_input('Read 2 fastq Location (R2.fastq): ')

    readList[read1] = read2

elif numFiles == 'n':
    others = 'Y'
    while others == 'Y':
        temp = raw_input('Read 1 fastq Location (R1.fastq/n if done): ')
        if not temp == 'n':
            read1 = temp
            read2 = raw_input('Read 2 fastq Location (R2.fastq): ')
            readList[read1] = read2
        elif temp == 'n':
            others = 'n'
        #others = raw_input('Enter more fastq files? (Y/n): ')

####################
#Load Previous Data#
####################
#previousDict = raw_input('Would you like to load previously sorted data? (Y/n): ')
if args.prevdict:
    prevDict = 'Y'
    prevDictLoc = args.prevdict
else:
    previousDict = 'n'
    prevDictLoc = 'null'

#previousDict = 'n'
#if previousDict == 'Y':
#    prevDictLoc = raw_input('Location of previous sorted data (/dir/data.pkl): ')
#elif previousDict == 'n':
#    prevDictLoc = 'null'

if not path.exists(outputDir):
    #make the output directory expanduser is used to allow ~/Desktop shortcuts
    mkdir(os.path.expanduser(outputDir))

####################
#Set Run Parameters#
####################
# set allowed umi mismatches
if args.umimismatch:
    distance_stringency = args.umimismatch
else:
    distance_stringency = 1

# variant threshold
if args.varthresh:
    varThresh = args.varthresh
else:
    varThresh = 0.75

# num supporting reads
if args.readsupport:
    supportingReads = args.readsupport
else:
    supportingReads = 5

# depth filtering
if args.dpfilter:
    DPNum = args.dpfilter
else:
    DPNum = 500

# ao filtering
if args.filterao:
    AONum = args.filterao
else:
    AONum = 5

# freebayes
if args.freebayes:
    freebayes = args.freebayes
else:
    freebayes = '/vol3/home/liggettl/programs/freebayes/bin/freebayes'

# read length
if args.readLength:
    readLength = args.readLength
else:
    readLength = False

#useDefaults = raw_input('Use Default Parameters? (Y/n): ')
#useDefaults = 'Y'
'''
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
'''
#create output directory
outputDir = outputDir + '/' + today + '_' + str(supportingReads) + '_' + str(varThresh)

# run locally or submit to cluster
if args.clustersubmit:
    clusterRun = 'Y'
else:
    clusterRun = 'n'
#clusterRun = raw_input('Submit to LRS cluster? (Y/n): ')
#clusterRun = 'n'

#estimate and output an error rate
if args.errorrate:
    errorRate = 'Y'
else:
    errorRate = 'n'

# trigger replacement of bad bases rather than elimination of entire capture
if args.badBaseSubstitute:
    badBaseSubstitute = True
else:
    badBaseSubstitute = False

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

#############################
#Record Files and Parameters#
#############################
def recordParams(parametersUsed, inputDir, read1, read2, varThresh, supportingReads, infoFile, REF, badBaseSubstitute):
    if infoOutput != 'n':
        info = open(infoFile, 'w')
        info.write(infoOutput)
        info.close()

    target = open(parametersUsed, 'w')
    target.write('Input File Location: %s\n' %(inputDir))
    target.write('HG Ref File Location: %s\n' %(REF))
    target.write("Read 1: %s\n" %(read1))
    target.write("Read 2: %s\n" %(read2))
    target.write("Distance Stringency: %d\n" %(distance_stringency))
    target.write("Variant Threshold: %f\n" %(varThresh))
    target.write("Supporting Reads: %d\n" %(supportingReads))
    target.write('DP Filter: %i\n' %(DPNum))
    target.write('AO Filter: %i\n' %(AONum))
    target.write('Read Length: %sbp\n' %(str(readLength)))
    if badBaseSubstitute:
        target.write('Bad Base Substitution Used\n')
    else:
        target.write('Bad Base Substitution Not Used\n')

    target.close()

############################
#Write Variables for Python#
############################
#outputs variables as pickle file that need to be loaded by cluster submitted
#scripts
def writePickle(one, two, specificOut):
    vardb = {}

    twoUmiOut = specificOut + '/twoUMIs.fastq'
    final_output_file = specificOut + '/finalOutput.fastq'
    coverage_file = specificOut + '/coverageData.txt'
    infoFile = specificOut + '/runInfo.txt'
    parametersUsed = specificOut + '/parametersUsed.txt'
    pickleOutput = specificOut + '/sortedSeqData.pkl'

    # this also records the parameters used in the run
    recordParams(parametersUsed, inputDir, read1, read2, varThresh, supportingReads, infoFile, REF, badBaseSubstitute)

    vardb['outputDir'] = specificOut
    vardb['varThresh'] = varThresh
    vardb['final_output_file'] = final_output_file
    vardb['supportingReads'] = supportingReads
    vardb['twoUmiOut'] = twoUmiOut
    vardb['distance_stringency'] = distance_stringency
    vardb['coverage_file'] = coverage_file
    vardb['previousDict'] = previousDict
    vardb['prevDictLoc'] = prevDictLoc
    vardb['pickleOutput'] = pickleOutput
    vardb['numFiles'] = numFiles
    vardb['read1'] = one
    vardb['read2'] = two
    vardb['inputDir'] = inputDir
    vardb['clusterRun'] = clusterRun
    vardb['noBigFiles'] = noBigFiles
    vardb['alignAndVar'] = alignAndVar
    vardb['dpnum'] = DPNum
    vardb['aonum'] = AONum
    vardb['freebayes'] = freebayes
    vardb['errorRate'] = errorRate
    vardb['readLength'] = readLength
    vardb['badBaseSubstitute'] = badBaseSubstitute
    vardb['REF'] = REF

    pickleVars = './variables.pkl'

    pickleFile = open(pickleVars, 'wb')
    pickle.dump(vardb, pickleFile)
    pickleFile.close()


#############
#Run main.py#
#############
#this requires that a shell script be created so that it can be sumitted
#to the LRS cluster not the python file
#look in ~/testScripts on Tesla for example
if __name__ == "__main__":

    # wait for any existing queueFile from still submitting jobs
    while os.path.exists('queueFile'):
        time.sleep(10)

    for i in readList:
        read1 = i
        read2 = readList[i]
        # make file specific output
        specificOut = outputDir + '/' + read1
        mkdir(os.path.expanduser(specificOut))

        while os.path.exists('queueFile'): # wait until previous file is processed before continuing
            time.sleep(10)
        writePickle(read1, read2, specificOut) # write new pickle with new read1/2
        system("touch queueFile") # set waiting file

        if clusterRun == 'n':
            system("python main.py")

        elif clusterRun == 'Y':
            system('bsub -n 1 < clusterSubmit.sh')
