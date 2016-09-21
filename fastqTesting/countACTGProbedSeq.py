#!/usr/bin/env python
'''
This script just counts the number of ATGC bases within a text
file. This was useful to count the number of bases in the
total probed sequence.
'''
import argparse
parser = argparse.ArgumentParser()

#####################
#Argparse Flag Input#
#####################
parser.add_argument('--indir', '-i', required=True, type=str, help='Specifies the input directory that contains the fastq files to be analyzed.')
parser.add_argument('--files', '-f', nargs='+', required=True, type=str, help='List of all fastq files to be analyzed.')

args = parser.parse_args()

##########
#Analysis#
##########
indir = args.indir
files = args.files

for i in files:
    bases = {'a':0, 't':0, 'g':0, 'c':0}
    inputFile = open(indir + i, 'r')
    for line in inputFile:
        for base in line:
            if base == 'a':
                bases['a'] += 1
            elif base == 't':
                bases['t'] += 1
            elif base == 'g':
                bases['g'] += 1
            elif base == 'c':
                bases['c'] += 1
    print('In file %s the base counts are:' % (i))
    print('A: %i\tT: %i\tG: %i\tC: %i\n' % (bases['a'],bases['t'],bases['g'],bases['c']))
