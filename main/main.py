import pickle
from os import system
from os import makedirs
from concatenateUMI import concatenateUMI
from goodCollapseDictionary import buildNestedDict
from goodCollapseDictionary import buildListDict
from goodCollapseDictionary import collapseNestedDict
from goodCollapseDictionary import collapseReadsListDict
from outputCoverage import outputCov
from time import time

start_time = time()

##################
#Unpack Variables#
##################
inputData = open('./variables.pkl', 'rb')
vardb = pickle.load(inputData)
inputData.close()
system('rm ./variables.pkl')
system('rm ./queueFile') # tell the system to start processing next file

read1 = vardb['read1']
read2 = vardb['read2']
outputDir = vardb['outputDir']
varThresh = vardb['varThresh']
final_output_file = vardb['final_output_file']
supportingReads = vardb['supportingReads']
twoUmiOut = vardb['twoUmiOut']
distance_stringency = vardb['distance_stringency']
coverage_file = vardb['coverage_file']
previousDict = vardb['previousDict']
prevDictLoc = vardb['prevDictLoc']
pickleOutput = vardb['pickleOutput']
numFiles = vardb['numFiles']
inputDir = vardb['inputDir']
clusterRun = vardb['clusterRun']
noBigFiles = vardb['noBigFiles']
alignAndVar = vardb['alignAndVar']
DPNum = vardb['dpnum']
AONum = vardb['aonum']
freebayes = vardb['freebayes']
errorRate = vardb['errorRate']

read1 = inputDir + '/' + read1
read2 = inputDir + '/' + read2

#####################
#Concatate R1 and R2#
#####################
concatenateUMI(read1, read2, twoUmiOut)

##############################
#Build/Get Seq Data Structure#
##############################
#build dict binning reads by concatenated UMIs
if previousDict == 'n':
    #seqDict = buildNestedDict(twoUmiOut, distance_stringency, pickleOutput)
    seqDict = buildListDict(twoUmiOut, distance_stringency, pickleOutput)

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
#collapseNestedDict(seqDict, varThresh, final_output_file, supportingReads, readLength)
averageErrorRate = collapseReadsListDict(seqDict, varThresh, final_output_file, supportingReads, readLength, errorRate)

#####################
#Output Seq Coverage#
#####################
outputCov(twoUmiOut, final_output_file, distance_stringency, coverage_file, averageErrorRate)

####################
#Align to Reference#
####################
if alignAndVar == 'Y':
    from align import align
    REF = '/vol3/home/liggettl/refgenomes/hg19.fa'
    bamOut = align(final_output_file, REF) # align and index

###############
#Call Variants#
###############
if alignAndVar == 'Y':
    from callVar import callVar
    vcfOut = bamOut.strip('bam') + 'vcf'
    callVar(freebayes, REF, bamOut, vcfOut)

###################
#VCF Decomposition#
###################
if alignAndVar == 'Y':
    from decomposeVCF import decompose
    blockDecomposedOut = decompose(vcfOut)

#################
#Filter Variants#
#################
# filters final vcf file to output either only AF=0 reads
# or AF=0.5 and AF=1 reads
if alignAndVar == 'Y':
    from varDPFilter import vcfFilter
    vcfFilter(inputDir, outputDir, blockDecomposedOut, AONum, DPNum)

################
#Output Runtime#
################
target = open(outputDir + '/runTime.txt', 'w')
target.write("Total Runtime:\n%s seconds" % (time() - start_time))
print("Total Runtime:\n%s seconds" % (time() - start_time))

##########################
#Remove Large Fastq Files#
##########################
# purpose of this is to suppress the output of large fastq files
# unless otherwise needed to prevent the use of unnecessary space
if noBigFiles == 'Y':
    system('rm %s' % (twoUmiOut))
