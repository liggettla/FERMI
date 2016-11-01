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
sample = inputDir + '/' + inSample + '/' + 'total_filtered.vcf'
#outDir = outputDir + '/' + inSample + '/'
outDir = outputDir + '/'

target = open(sample, 'r')

#################
# Datastructure #
#################
from copy import deepcopy
#varDict = {'c': {'a':1, 't':2 ...} ...}
baseDict = {'A':0, 'T':0, 'G':0, 'C':0}
# if copy() is not used, only one instance of the baseDict is used
# under each of the keys, and cannot be independently changed
varDict = {'A':baseDict.copy(), 'T':baseDict.copy(), 'G':baseDict.copy(), 'C':baseDict.copy()}
fullDict = {'oncoSites':deepcopy(varDict), 'oncoGenes':deepcopy(varDict), 'TIIIRegions':deepcopy(varDict)}

###################
# Classify Sample #
###################
# determine if sample is oncogenic, non-oncogenic but in an oncogene or TIII
oncoSites = [5073770, 7577539,7577119,115256529,115258747,115258744,534287,534288,534289,25398284,25380275,106197266,106197267,106197268,106197269,106155172,106155173,106155174,25457242,25457243,209113112,209113113,90631934,90631838,48649700]
oncoGenes = []
for i in oncoSites: # provides total oncogene region (may underestimate total size)
    i = i / 1000
    oncoGenes.append(i)
TIIIRegions = ['115227', '229041', '110541', '112997', '121167', '123547', '124428', '1397', '2126', '2390', '2593', '11486', '92527', '73379', '82455', '85949']

###############
# Read Sample #
###############
a = 'A'
t = 'T'
c = 'C'
g = 'G'
baseList = [a, t, g, c]


for i in target:
    if '#' not in i and 'chr' in i: #skip the damn info
        ref = i.split('\t')[3]
        alt = i.split('\t')[4]
        loc = i.split('\t')[1]

        #import pdb
        #pdb.set_trace()
        for j in baseList:
            if j == ref:
                for k in baseList:
                    if k == alt:

                        if int(loc) in oncoSites: #is the var oncogenic?
                            site = 'oncoSites'
                            fullDict['oncoSites'][j][k] += 1
                        elif int(loc)/1000 in oncoGenes: #is the var nonOncogenic exomic?
                            site = 'oncoGenes'
                            fullDict['oncoGenes'][j][k] += 1
                        elif str(int(loc)/1000) in TIIIRegions: #is the var TIII?
                            site = 'TIIIRegions'
                            fullDict['TIIIRegions'][j][k] += 1

##################
# Export Results #
##################
from os import system
#system('mkdir %s' % (outDir))
for region in fullDict:
    output = outDir + region + '.txt' # output individual files for oncoGenes, TIII, oncoSites
    with open(output, 'w') as out:
        out.write('Ref\tVar\tObs\n') # headers
        for ref in fullDict[region]:
            for var in fullDict[region][ref]:
                num = fullDict[region][ref][var]
                out.write('%s\t%s\t%i\n' % (ref, var, num))
        out.close()
