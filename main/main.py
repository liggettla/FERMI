import pickle
from os import system
from os import makedirs
from concatenateUMI import concatenateUMI
from concatenateUMI import duplexConcatenate
from goodCollapseDictionary import buildListDict
from goodCollapseDictionary import collapseReadsListDict
from goodCollapseDictionary import duplexCollapse
from goodCollapseDictionary import get_one_bp_mismatches
from goodCollapseDictionary import find_complementary_umis
from goodCollapseDictionary import collapse_paired_reads
from goodCollapseDictionary import calcCoverageError
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
readLength = vardb['readLength']
badBaseSubstitute = vardb['badBaseSubstitute']
REF = vardb['REF']
duplex = vardb['duplex']
minimalOutput = vardb['minimalOutput']
realvsmock = vardb['realvsmock']

read1 = inputDir + '/' + read1
read2 = inputDir + '/' + read2

#####################
#Concatate R1 and R2#
#####################
if not duplex:
    print('Concatenating UMIs...')
    concatenateUMI(read1, read2, twoUmiOut)
elif duplex:
    print('Duplex Concatenating...')
    duplexConcatenate(read1, read2, twoUmiOut, realvsmock)

##############################
#Build/Get Seq Data Structure#
##############################
#build dict binning reads by concatenated UMIs
print('Building Dictionary...')
if previousDict == 'n':
    seqDict = buildListDict(twoUmiOut, distance_stringency, pickleOutput)

#retrieve previous seq data structure
elif previousDict == 'Y':
    prevData = open(prevDictLoc, 'rb')
    seqDict = pickle.load(prevData)
    prevData.close()

####################
# Calc Read Length #
####################
#calculate individual read length between the UMIs
print('Computing Read Lengths...')
with open(twoUmiOut, 'r') as target:
    header = next(target)
    readSeq = next(target).rstrip('\n')

    # if readLength not set use full seq, else use specified num
    if not readLength:
        readLength = len(readSeq) - 12

################
#Collapse Reads#
################
#collapse reads on binned UMI data structure
print('Collapsing Reads...')
averageErrorRate, averageCoverage = collapseReadsListDict(seqDict, varThresh, final_output_file, supportingReads, readLength, errorRate, badBaseSubstitute)

###################
# Duplex Collapse #
###################
# Duplex collapse using the two initial two strands of every capture
# to eliminate any dissimilar variants
'''
print('Duplex Collapsing Reads...')
if duplex:
    #averageErrorRate, averageCoverage = duplexCollapse(seqDict, varThresh, final_output_file, supportingReads, readLength, errorRate, badBaseSubstitute)
    coverageList, errorRateList, duplexDict = duplexCollapse(seqDict, varThresh, final_output_file, supportingReads, readLength, errorRate, badBaseSubstitute)
    deDuplexList = find_complementary_umis(duplexDict)
    collapse_paired_reads(deDuplexList, readLength, final_output_file)
    averageErrorRate, averageCoverage = calcCoverageError(coverageList, errorRateList, errorRate)
'''

#####################
#Output Seq Coverage#
#####################
print('Computing Coverages...')
print('Herding Llamas...')
outputCov(twoUmiOut, final_output_file, distance_stringency, coverage_file, averageErrorRate, averageCoverage)

####################
#Align to Reference#
####################
if alignAndVar == 'Y':
    print('Aligning to Reference...')
    from align import align
    bamOut = align(final_output_file, REF) # align and index

###############
#Call Variants#
###############
if alignAndVar == 'Y':
    print('Calling Variants...')
    from callVar import callVar
    vcfOut = bamOut.strip('bam') + 'vcf'
    callVar(freebayes, REF, bamOut, vcfOut)

###################
#VCF Decomposition#
###################
if alignAndVar == 'Y':
    print('Decomposing VCFs...')
    from decomposeVCF import decompose
    blockDecomposedOut = decompose(vcfOut)

#################
#Filter Variants#
#################
# filters final vcf file to output either only AF=0 reads
# or AF=0.5 and AF=1 reads
if alignAndVar == 'Y':
    print('Filtering Allele Frequencies...')
    from varDPFilter import vcfFilter
    vcfFilter(inputDir, outputDir, blockDecomposedOut, AONum, DPNum)

##########################
#Remove Large Fastq Files#
##########################
# purpose of this is to suppress the output of large fastq files
# unless otherwise needed to prevent the use of unnecessary space
if noBigFiles == 'Y':
    print('Cleaning up Fastqs...')
    system('rm %s' % (twoUmiOut))

##############################
# Remove Non-Probed Variants #
##############################
# the purpose of this is to create a final vcf file that only contains
# those variants that fall within probed regions of the genome
if alignAndVar == 'Y':
    print('Eliminating Off-Target Alignments...')
    from eliminateNonspecificReads import elimBadAligns
    unFiltered = outputDir + '/total_filtered.vcf'
    filtered = outputDir + '/onlyProbedRegions.vcf'
    elimBadAligns(unFiltered, filtered)

############################
# Plot Mutations Per Probe #
############################
# the purpose of this is to understand if mutations cluster within
# particular probes more than others
if alignAndVar == 'Y':
    print('Computing Probe Bias...')
    from mutationsPerProbe import mutationsPerProbe
    mutationsPerProbe(filtered, outputDir)

##########################
# Remove Undesired Files #
##########################
# this will remove many of the output files that are mostly
# used for troubleshooting and not for typical analysis
if minimalOutput:
    print('Cleaning up Outputs...')
    system('rm %s' % (outputDir + '/finalOutput.fastq'))
    system('rm %s' % (outputDir + '/finalOutput.bam'))
    system('rm %s' % (outputDir + '/finalOutputBlockDecomposed.vcf'))
    system('rm %s' % (outputDir + '/finalOutputDecomposed.vcf'))
    system('rm %s' % (outputDir + '/finalOutput.vcf'))
    system('rm %s' % (outputDir + '/total_filtered.vcf'))
    system('rm %s' % (outputDir + '/AF0_filtered.vcf'))
    system('rm %s' % (outputDir + '/AF1_filtered.vcf'))
    system('rm %s' % (outputDir + '/finalOutput.bam'))

################
#Output Runtime#
################
# keep this as the final task that is run
print('Computing Runtime...')
target = open(outputDir + '/runTime.txt', 'w')
target.write("Total Runtime:\n%s seconds" % (time() - start_time))
print("Total Runtime:\n%s seconds" % (time() - start_time))
print('Analysis Complete.\n')
