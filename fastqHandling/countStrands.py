#!/usr/bin/env python
# The purpose of this script is to count within original fastq
# files how many occurrances there are of particular variants
# in order to look for linked variants

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', type=str, help='Input file directory containing original fastq files.')

args = parser.parse_args()
indir = args.input

import os

def getFileList():
    for file in os.listdir(indir):
        if file.endswith('.fastq'):
            inFile = indir + '/' + file
            countVars(inFile)

def countVars(inFile):
    target = open(inFile, 'r')
    possibleSeqs = {0:'ATACCCAAGGGAGAGTCCCTGCTACAATACAGGTTGTGGCATTATAGAAGACTAAAGTAGGAGTGACATAGAAGATAAGTATATCCAAGTATATCCATAG',
           1:'ATACCCAAGGGAGAGTCCCTGCTACAATACAGGTTGTGGCATTATAGAAGACTAAAGTAGGAGCGACATAGAAGATAAGTATATCCAAGTATATCCATAG',
           2:'ATACCCAAGGGAGAGTCCCTGCTACAATACAGGTTGTGGCATTATAGAAGACTAAAGTAGGAGTGACACAGAAGATAAGTATATCCAAGTATATCCATAG',
           3:'ATACCCAAGGGAGAGTCCCTGCTACAATACAGGTTGTGGCATTATAGAAGACTAAAGTAGGAGTGACATAGAAGATAAGTATATCCAAGTATATCCACAG',
           12:'ATACCCAAGGGAGAGTCCCTGCTACAATACAGGTTGTGGCATTATAGAAGACTAAAGTAGGAGCGACACAGAAGATAAGTATATCCAAGTATATCCATAG',
           13:'ATACCCAAGGGAGAGTCCCTGCTACAATACAGGTTGTGGCATTATAGAAGACTAAAGTAGGAGCGACATAGAAGATAAGTATATCCAAGTATATCCACAG',
           23:'ATACCCAAGGGAGAGTCCCTGCTACAATACAGGTTGTGGCATTATAGAAGACTAAAGTAGGAGTGACACAGAAGATAAGTATATCCAAGTATATCCACAG',
           123:'ATACCCAAGGGAGAGTCCCTGCTACAATACAGGTTGTGGCATTATAGAAGACTAAAGTAGGAGCGACACAGAAGATAAGTATATCCAAGTATATCCACAG'}
    counts = {0:0, 1:0, 2:0, 3:0, 12:0, 13:0, 23:0, 123:0}

    for line in target:
        for seq in possibleSeqs:
            if possibleSeqs[seq] in line:
                counts[seq] += 1

    for i in counts:
        print str(i) + ':' + str(counts[i])

if __name__ == "__main__":
    getFileList()
