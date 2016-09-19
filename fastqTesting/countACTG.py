#!/usr/bin/env python
'''
This script tests the ACTG content of the fastq file.
This script was a suggestion to try and figure out why there
is an unmatched C-T G-C frequency.
'''
import argparse
parser = argparse.ArgumentParser()

#####################
#Argparse Flag Input#
#####################
parser.add_argument('--indir', '-i', required=True, type=str, help='Specifies the input directory that contains the fastq files to be analyzed.')
parser.add_argument('--outdir', '-o', required=True, type=str, help='Point to the output directory')
parser.add_argument('--files', '-f', nargs='+', required=True, type=str, help='List of all fastq files to be analyzed.')

args = parser.parse_args()

##########
#Analysis#
##########
indir = args.indir
outdir = args.outdir
fastqList = args.files
outputFile = open(outdir + 'basecounts', 'w')

for i in fastqList:
    bases = {'a':0, 't':0, 'g':0, 'c':0}
    inputFile = open(indir + i, 'r')

    position = 1
    for line in inputFile:
	if position == 1:
	    position += 1
	elif position == 2:
	    readSeq = line[6:-6]
            for base in readSeq:
                if base == 'A':
                    bases['a'] += 1
                elif base == 'T':
                    bases['t'] += 1
                elif base == 'G':
                    bases['g'] += 1
                elif base == 'C':
                    bases['c'] += 1
	    position += 1
	elif position == 3:
	    position += 1
	elif position == 4:
	    position = 1

    outputFile.write('Input File:\n%s\n' % (indir + i))
    outputFile.write('A: %i\tT: %i\tG: %i\tC: %i\n\n' % (bases['a'], bases['t'], bases['g'], bases['c']))

