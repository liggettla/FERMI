# the purpose of this script is to compare two repeatedly sequenced samples originating
# from the same DNA sample to understand how repeatable the VAF is for a given mutation
# in each of the analysis files
# for example in 1.vcf VAF of a variant might be 0.5 but in 2.vcf it is 0.2

############
# Argparse #
############
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--indir', '-i', type=str, required=True, help='Input directory containing the vcf files to be analyzed')
parser.add_argument('--outdir', '-o', type=str, required=True, help='Output directory for plots')
parser.add_argument('--samples', '-s', type=str, nargs='*', required=True, help='Filenames of the two samples to be compared')

args = parser.parse_args()

inputDir = args.indir
outputDir = args.outdir
samples = args.samples # this is a list

#################
# Read In Files #
#################
from itertools import combinations
sample1 = inputDir + '/' + samples[0]
sample2 = inputDir + '/' + samples[1]

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
            loc = str(i.split('\t')[1])
            AO = i.split(';')[5]
            DP = i.split(';')[7]
            var = i.split('\t')[4] # the observed bp change

            AONum = float(AO.split(',')[0][3:])
            DPNum = float(DP.split(',')[0][3:])
            AFNum = AONum / DPNum

            dataframe[loc] = {'var': var, 'vaf': AFNum}

###################
# Get Common Vars #
###################
outputFile = outputDir + '/vafRepeatability.txt'
output = open('outputFile', 'w')
output.write('Sample1\tSample2\n')

for i in df1:
    if i in df2:
        output.write('%s\t%s\n' % (df1[i]['vaf'], df2[i]['vaf']))

output.close()

################
# Move Results #
################
from os import system
command = 'mv outputFile %s' % (outputFile)
system(command)





