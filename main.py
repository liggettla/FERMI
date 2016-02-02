import pickle
from os import system
from concatenateUMI import concatenateUMI
from concatenateUMI import concatenateUMI
from goodCollapseDictionary import buildNestedDict
from outputCoverage import outputCov
from goodCollapseDictionary import buildNestedDict
from goodCollapseDictionary import collapseNestedDict

##################
#Unpack Variables#
##################
inputData = open('./variables.pkl', 'rb')
vardb = pickle.load(inputData)
inputData.close()
system('rm ./variables.pkl')

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

#####################
#Align and Call Vars#
#####################
REF = '/vol3/home/liggettl/refgenomes/hg19.fa'
bamOut = final_output_file.strip('fastq') + 'bam'

system("bwa mem %s %s | samtools view -bS - | samtools sort > %s" % (REF, final_output_file, bamOut))
system("samtools index %s" % (bamOut))

'''
Nothing below this has yet been implemented
#######################
#Output Var Freq Table#
#######################
outputDir = outputDir + '/'
vcfFile = outputDir + '/finalOutput.vcf'
system('cp %s ./' % (vcfFile))
system('bash newIdentifyVars.sh')
system('mv finalOutput.vcf %s' (outputDir))

###################
#Plot Allele Freqs#
###################
system('Rscript plotVarFreq.R')
system('mv allelefreqs.txt %s' (outputDir))
system('mv allelefreqs.jpg %s' (outputDir))


#this should output a plot automatically, but is not runtime tested
#does not output log10() of plot yet
#it may be best to just system(Rscript plotting.R) somehow
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
'''
