#!/usr/bin/env python
# The purpose of this script is to quantify the base changes by base
# in order to understand if there is a bias to a particular base change
# as the observed C to T change

############
# Argparse #
############
import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--indir', '-i', required=True, type=str, help='Specifies the input directory containing all folders containing output analysis from a fermi analysis run.')
parser.add_argument('--outdir', '-o', required=True, type=str, help='Specifies the output directory location for the analysis output file')
parser.add_argument('--sample', '-s', required=True, type=str, help='Name of the directory containing fermi analysis of sample to be compared.')

args = parser.parse_args()

inputDir = args.indir
outputDir = args.outdir
inSample = args.sample
sample = inputDir + '/' + inSample + '/' + 'finalOutputBlockDecomposed.vcf'

target = open(sample, 'r')

#################
# Datastructure #
#################
#varDict = {'c': {'a':1, 't':2 ...} ...}
baseDict = {'A':0, 'T':0, 'G':0, 'C':0}
# if copy() is not used, only one instance of the baseDict is used
# under each of the keys, and cannot be independently changed
varDict = {'A':baseDict.copy(), 'T':baseDict.copy(), 'G':baseDict.copy(), 'C':baseDict.copy()}

###############
# Read Sample #
###############
a = 'A'
t = 'T'
c = 'C'
g = 'G'
baseList = [a, t, g, c]

#import pdb
#pdb.set_trace()

for i in target:
    if '#' not in i and 'chr' in i: #skip the damn info
        ref = i.split('\t')[3]
        alt = i.split('\t')[4]

        for j in baseList:
            if j == ref:
                for k in baseList:
                    if k == alt:
                        varDict[j][k] += 1
from pprint import pprint
pprint(varDict)

