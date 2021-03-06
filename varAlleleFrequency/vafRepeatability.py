#!/usr/bin/python
# the purpose of this script is to compare two repeatedly sequenced samples originating
# from the same DNA sample to understand how repeatable the VAF is for a given mutation
# in each of the analysis files
# for example in 1.vcf VAF of a variant might be 0.5 but in 2.vcf it is 0.2

############
# Argparse #
############
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--indir', '-i', type=str, required=True, help='Input directory containing the vcf files to be analyzed: /dir')
parser.add_argument('--outdir', '-o', type=str, required=True, help='Output directory for plots: /dir')
parser.add_argument('--samples', '-s', type=str, nargs='*', required=True, help='Filenames of the two samples to be compared: A1-R1.fastq')
#parser.add_argument('--rarevars', '-r', action='store_true', help='Only include rare variants in the analysis that are at an AF<0.30')
parser.add_argument('--rarevars', '-r', type=float, help='This can be set to cutoff the data at a certain allele frequency and only include variants below a particular frequency like 0.03 or 0.003.')
parser.add_argument('--commonVars', '-c', action='store_true', help='This will only plot variants that are found in both samples and ignore those variants that are only found in one of the samples.')
parser.add_argument('--germline', '-g', type=str, nargs='*', help='Only output those variants that changed from these bases.')
parser.add_argument('--variant', '-v', type=str, nargs='*', help='Only output those variants that change to these bases.')

args = parser.parse_args()

inputDir = args.indir
outputDir = args.outdir
samples = args.samples # this is a list

if args.germline:
    germline = args.germline
else:
    germline = ['A','T','G','C']

if args.variant:
    variant = args.variant
else:
    variant = ['A','T','G','C']

#################
# Read In Files #
#################
from itertools import combinations
sample1 = inputDir + '/' + samples[0] + '/' + 'total_filtered.vcf'
sample2 = inputDir + '/' + samples[1] + '/' + 'total_filtered.vcf'
#sample1 = inputDir + '/' + samples[0] + '/' + 'onlyProbedRegions.vcf'
#sample2 = inputDir + '/' + samples[1] + '/' + 'onlyProbedRegions.vcf'

file1 = open(sample1, 'r')
df1 = {}
file2 = open(sample2, 'r')
df2 = {}

# iterate through input files
for x in range(1,3):
    if x == 1:
        target = file1
        dataframe = df1
    if x == 2:
        target = file2
        dataframe = df2

# populate dataframes
    for i in target:
        if '#' not in i and 'chr' in i: # skip the info
            chrom = str(i.split('\t')[0])
            loc = str(i.split('\t')[1])
            AO = i.split(';')[5]
            DP = i.split(';')[7]
            var = i.split('\t')[4] # the observed bp change
            WT = i.split('\t')[3] # WT base

            AONum = float(AO.split(',')[0][3:])
            DPNum = float(DP.split(',')[0][3:])
            AFNum = AONum / DPNum

            if args.rarevars:
                cutoff = args.rarevars
                if AFNum < cutoff:
                    dataframe[loc] = {'var': var, 'vaf': AFNum, 'chr': chrom, 'wt': WT}
            else:
                dataframe[loc] = {'var': var, 'vaf': AFNum, 'chr': chrom, 'wt': WT}

###################
# Get Common Vars #
###################
outputFile = outputDir + '/vafRepeatability.txt'
plotFile1 = outputDir + '/vafRepeatabilityRegression.jpg'
plotFile2 = outputDir + '/vafRepeatabilityNoRegression.jpg'
output = open('outputFile', 'w')
output.write('Sample1\tSample2\tIdentity\n')

# only output variants observed in both samples
if args.commonVars:
    for i in df1:
        if i in df2:
            if df1[i]['wt'] in germline and df1[i]['var'] in variant:
                output.write('%s\t%s\t%s:%s:%s->%s\n' % (df1[i]['vaf'], df2[i]['vaf'], df1[i]['chr'], i, df1[i]['wt'], df1[i]['var']))

# output all variants including those not observed in both samples
else:
    for i in df1:
        if df1[i]['wt'] in germline and df1[i]['var'] in variant:
            # write overlapping variants
            if i in df2:
                output.write('%s\t%s\t%s:%s:%s->%s\n' % (df1[i]['vaf'], df2[i]['vaf'], df1[i]['chr'], i, df1[i]['wt'], df1[i]['var']))
            # write variants found only in df1
            else:
                output.write('%s\t%s\t%s:%s:%s->%s\n' % (df1[i]['vaf'], 0, df1[i]['chr'], i, df1[i]['wt'], df1[i]['var']))

    # write variants found only in df2
    for i in df2:
        if df2[i]['wt'] in germline and df2[i]['var'] in variant:
            if i not in df1:
                output.write('%s\t%s\t%s:%s:%s->%s\n' % (0, df2[i]['vaf'], df2[i]['chr'], i, df2[i]['wt'], df2[i]['var']))

output.close()

#############
# Plot Data #
#############
from os import system
command = 'Rscript plotvafRepeatability.R'
system(command)

################
# Move Results #
################
command = 'mv outputFile %s' % (outputFile)
system(command)
command = 'mv output1.jpg %s' % (plotFile1)
system(command)
command = 'mv output2.jpg %s' % (plotFile2)
system(command)

#############
# Plot Data #
#############
'''
command = 'eog vafRepeatabilityRegression.jpg'
system(command)
command = 'eog vafRepeatabilityNoRegression.jpg'
system(command)
'''
command = 'eog %s' % (plotFile2)
system(command)


