#!/usr/bin/python

# This script provides an algorithm for automated sorting
# through FERMI pipeline analyzed files for analysis.
# This will be useful for automated selection of single
# samples to compare against grouped others.

############
# Argparse #
############
import argparse
from numpy import mean

parser = argparse.ArgumentParser()
parser.add_argument('--indir', '-i', type=str, required=True, help='Input directory containing the vcf files to be analyzed: /dir')
parser.add_argument('--outdir', '-o', type=str, required=True, help='Output directory for plots: /dir')
parser.add_argument('--principle', '-p', type=str, required=True, help='This is the principle sample being compared to an averaged set of other samples. Ex: A1-R1')


from os import listdir
fileList = []
for i in listdir(inDir):
    fileList.append(i)

for i in fileList:
    print i
