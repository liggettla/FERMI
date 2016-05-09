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

###################
#Read Through VCFs#
###################
sample1 = inputDir + '/' + samples[0] + '/' + 'finalOutputBlockDecomposed.vcf'
sample2 = inputDir + '/' + samples[0] + '/' + 'finalOutputBlockDecomposed.vcf'
file1 = open(sample1, 'r')
file2 = open(sample2, 'r')

for i in file1:
    presence = 'n'

    lowFreqPresent = 0
    lowFreqAbsent = 0
    highFreqPresent = 0
    highFreqAbsent = 0

    if '#' not in i and 'chr' in i: #skip the damn info

        loc = i.split('\t')[1]
        AO = i.split(';')[5]
        DP = i.split(';')[7]
        var = i.split('\t')[4]

        AONum = float(AO.split(',')[0][3:])
        DPNum = float(DP.split(',')[0][3:])
        AFNum = AONum / DPNum

        for j in file2:
            loc2 = i.split('\t')[1]
            var2 = i.split('\t')[4]

            if loc == loc2 and var == var2:
                presence = 'y'

        if presence == 'y':
            if AFNum < 0.5:
                lowFreqPresent += 1
            else:
                highFreqPresent += 1
        else:
            if AFNum < 0.5:
                lowFreqAbsent += 1
            else:
                highFreqAbsent += 1


    outTarget.write('Low Freq Var Present: \n' + str(lowFreqPresent) + '\n')
    outTarget.write('High Freq Var Present: \n' + str(highFreqPresent) + '\n')
    percentLow = float(lowFreqPresent) / (lowFreqAbsent + lowFreqPresent)
    outTarget.write('Percent Low Freq Var Present: \n' + percentLow + '\n')
    percentHigh = float(highFreqPresent) / (highFreqAbsent + highFreqPresent)
    outTarget.write('Percent High Freq Var Present: \n' + percentHigh + '\n')





