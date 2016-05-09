#!/usr/bin/env python

# The purpose of this script is to understand how well mutations are
# repeatedly seen between multiple sequence runs of DNA originating from
# the same source.
import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--indir', '-i', required=True, type=str, help='Specifies the input directory containing all folders containing output analysis from a fermi analysis run.')
parser.add_argument('--outdir', '-o', required=True, type=str, help='Specifies the output directory location for the analysis output file')
parser.add_argument('--samples', '-s', required=True, nargs='*', type=str, help='Name of the directories containing fermi analysis of samples to be compared.')

args = parser.parse_args()

inputDir = args.indir
outputDir = args.outdir
samples = args.samples # this is a list

outputFile = outputDir + '/' + 'sampleRepeatability.txt'
outTarget = open(outputFile, 'w')

###############
#Variant Class#
###############
class Variant:
    def __init__(self, location, AO, DP, variant, AF):
        self.location = location
        self.AO = int(AO)
        self.DP = int(DP)
        self.variant = variant
        self.AF = AF

#################
#Build Dataframe#
#################
# {sample1:[var1, var2, var3], sample2: [var1, var2, var3]}
sampleDict={}

for sample in samples:
    varList = []
    inFile = inputDir + '/' + sample + '/' + 'finalOutputBlockDecomposed.vcf'
    target = open(inFile, 'r')
    for line in target:
        if '#' not in line and 'chr' in line: #skip the damn info
            loc = line.split('\t')[1]
            AO = line.split(';')[5]
            DP = line.split(';')[7]
            var = line.split('\t')[4]

            AONum = float(AO.split(',')[0][3:])
            DPNum = float(DP.split(',')[0][3:])
            AFNum = AONum / DPNum

            variant = Variant(loc, AONum, DPNum, var, AFNum)
            varList.append(variant)

    sampleDict[str(sample)] = varList
    #print sampleDict[sample][24228].DP


    target.close()

###############
#Compare Files#
###############
from itertools import combinations

for a, b in combinations(samples, 2):
    sample1 = sampleDict[a]
    sample2 = sampleDict[b]
    len1 = len(sample1)
    len2 = len(sample2)

    lowFreqPresent = 0
    lowFreqAbsent = 0
    highFreqPresent = 0
    highFreqAbsent = 0


    for i in range(len1):
        presence = 'n'
        loc1 = sampleDict[a][i].location
        var1 = sampleDict[a][i].variant
        AF1 = sampleDict[a][i].AF
        AO1 = sampleDict[a][i].AO

        for j in range(len2):
            loc2 = sampleDict[b][i].location
            var2 = sampleDict[b][i].variant
            AF2 = sampleDict[b][i].AF
            AO2 = sampleDict[b][i].AO

            if loc1 == loc2 and var1 == var2:
                presence = 'y'

        if presence == 'y':
            if AF1 == 0:
                lowFreqPresent += 1
            else:
                highFreqPresent += 1
        else:
            if AF1 == 0:
                lowFreqAbsent += 1
            else:
                highFreqAbsent += 1

    outTarget.write('Low Freq Var Present: \n' + str(lowFreqPresent) + '\n')
    outTarget.write('High Freq Var Present: \n' + str(highFreqPresent) + '\n')
    percentLow = float(lowFreqPresent) / (lowFreqAbsent + lowFreqPresent)
    outTarget.write('Percent Low Freq Var Present: \n' + percentLow + '\n')
    percentHigh = float(highFreqPresent) / (highFreqAbsent + highFreqPresent)
    outTarget.write('Percent High Freq Var Present: \n' + percentHigh + '\n')














'''

######################
#Build Donor Loc List#
######################
donorList = []
target = open(donorFile, 'r')
for line in target:
    if '#' not in line and 'chr' in line: #skip the damn info
            loc = line.split('\t')[1]
            donorList.append(loc)
target.close()

##########################
#Build Recipient Loc Dict#
##########################
recAF0Dict = {}
recAF1Dict = {}
for i in recSamples:
    A1 = inputDir + '/' + i + '/AF1_filtered.vcf'
    target1 = open(A1, 'r')

    AF1List = []
    for line in target1:
        if '#' not in line and 'chr' in line: #skip the damn info
            loc = line.split('\t')[1]
            AF1List.append(loc)

    recAF1Dict[i] = AF1List

target1.close()

#####################
#Get Uniq Donor Vars#
#####################
commonDonor = []
uniqDonor = []
for i in donorList:
    for j in recAF1Dict:
        if i in recAF1Dict[j]:
            commonDonor.append(i)

for i in donorList:
    if not i in commonDonor:
        uniqDonor.append(i)

#####################
#Output Observations#
#####################
results = inputDir + '/dilutionResults.txt'
target = open(results, 'w')

for i in recSamples:
    A0 = inputDir + '/' + i + '/AF0_filtered.vcf'
    target0 = open(A0, 'r')

    for line in target0:
        if '#' not in line and 'chr' in line: #skip the damn info
            loc = line.split('\t')[1]
            AO = line.split(';')[5]
            DP = line.split(';')[7]
            var = line.split('\t')[4]
            #AF = line.split(';')[3]

            AONum = float(AO.split(',')[0][3:])
            DPNum = float(DP.split(',')[0][3:])
            AFNum = AONum / DPNum

            if loc in uniqDonor:
                target.write('Sample: ' + str(i) + '\n' + 'Location: ' + str(loc) + '\n' \
                        + 'AO: ' + str(AONum) + '\n' \
                        + 'DP: ' + str(DPNum) + '\n' \
                        + 'AF: ' + str(AFNum) + '\n')

target.close()
target0.close()





recipientCons = 'AF0_filtered.vcf'



consName = '_finalOutput.vcf'

spikeList = []
List24 = []
List25 = []
List28 = []
List29 = []
listOfLists = [List24, List25, List28, List29]

target = open(inputDir + donorSample + consName, 'r')

for line in target:
    if '#' not in line and 'chr' in line: #skip the damn info
        AO = line.split(';')[5]
        DP = line.split(';')[7]
        AF = line.split(';')[3]

        AONum = int(AO.split(',')[0][3:])
        AFNum = float(AF.split(',')[0][3:])
        DPNum = int(DP.split(',')[0][3:])

#get spike-in vars
        if AONum > 5 and (AFNum == 1 or AFNum == 0.5) and DPNum > 50:
            loc = line.split('\t')[1]
            spikeList.append(loc)

#get all variants
for i in listOfLists:
    index = listOfLists.index(i)
    target = open(inputDir + samples[index] + consName, 'r')
    for line in target:
        if '#' not in line and 'chr' in line: #skip the damn info
            loc = line.split('\t')[1]
            i.append(loc)

outFile = open(inputDir + 'DilutionEfficiency.txt', 'w')
numVars = len(spikeList)

for sample in listOfLists:
    counter = 0
    index = listOfLists.index(sample)
    for loc in sample:
        if loc in spikeList:
            counter += 1
    fraction = float(counter) / numVars * 100
    outFile.write('In sample %s %f percent of donor in variants were observed.\n' % (samples[index], fraction))

print len(spikeList)

'''
