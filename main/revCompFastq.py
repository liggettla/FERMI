'''
This script converts fastq files into their reverse
complements. The reason for doing this is as a quality control
to see if the particular frequencies of base changes are
reversed. If some systematic error exists that is causing
more C-T changes in the code, then this might also appear in
the reverse complement.
'''
from Bio.Seq import Seq
import argparse

############
# Argparse #
############
parser = argparse.ArgumentParser()
parser.add_argument('-indir', '-i', type=str, required=True, help='Specifies the input directory.')
parser.add_argument('-outdir', '-o', type=str, required=True, help='Specifies the output directory.')
parser.add_argument('-files', '-f', type='+', help='List of input files.')

args = parser.parse_args()

##################
# Rev Compliment #
##################
indir = args.indir
outdir = args.outdir
files = args.files
position = 1
for file in
inFile = open(,'r')
