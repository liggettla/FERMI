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

            outputFile = indir + '/' + 'tripleSNPCounts.txt'
            target = open(outputFile, 'a')

            countVars(inFile, file)

def countVars(inFile, file):
    outputFile = indir + '/' + 'tripleSNPCounts.txt'
    target = open(outputFile, 'a')
    target.write(file)
    print '\n'
    print file
    target = open(inFile, 'r')

    # These are vars that are possibly linked within the 1.12.2016 dataset
    possibleSeqs = {0:'TGACATAGAAGATAAGTATATCCAAGTATATCCATAG',
           1:'CGACATAGAAGATAAGTATATCCAAGTATATCCATAG',
           2:'TGACACAGAAGATAAGTATATCCAAGTATATCCATAG',
           3:'TGACATAGAAGATAAGTATATCCAAGTATATCCACAG',
           12:'CGACACAGAAGATAAGTATATCCAAGTATATCCATAG',
           13:'CGACATAGAAGATAAGTATATCCAAGTATATCCACAG',
           23:'TGACACAGAAGATAAGTATATCCAAGTATATCCACAG',
           123:'CGACACAGAAGATAAGTATATCCAAGTATATCCACAG'}
    counts = {0:0, 1:0, 2:0, 3:0, 12:0, 13:0, 23:0, 123:0}

    # Scan through fastq file and look for occurances of above seqs
    for line in target:
        for seq in possibleSeqs:
            if possibleSeqs[seq] in line:
                counts[seq] += 1

    for i in counts:
        target = open(outputFile, 'a')
        print str(i) + ':' + str(counts[i])
        target.write(str(i) + ':' + str(counts[i]) + '\n')
        target.close()

if __name__ == "__main__":
    getFileList()
