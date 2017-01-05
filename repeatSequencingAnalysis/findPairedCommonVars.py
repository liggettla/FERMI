#!/usr/bin/python

# The purpose of this script is to take an input vcf and display
# variants that are above certain thresholds to allow the user to
# identify homo/het variants that are found within the same UMI capture

from parseLine import seqRead

def runArgParse():
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--infile', '-i', type=str, required=True, help='The input vcf file to be analyzed, Ex: /file.vcf')
    parser.add_argument('--afcutoff', '-a', type=float, required=True, help='The inclusive allele frequency cutoff below which variants will not be displayed.')
    parser.add_argument('--aocutoff', '-o', type=int, required=True, help='The inclusive coverage cutoff below which variants will not be displayed.')

    args = parser.parse_args()

    inFile = args.infile
    afCutoff = args.afcutoff
    aoCutoff = args.aocutoff

    return inFile, afCutoff, aoCutoff

# Read through input vcf file
def readFile(inFile, afCutoff, aoCutoff):
    target = open(inFile, 'r')

    print 'Chrom-Loc\tAO\t\tDP\t\tAF'
    for line in target:
	if '#' not in line and 'chr' in line: # skip vcf info
            lineObj = seqRead(line)
            if (lineObj.af() >= afCutoff) and (lineObj.ao() >= aoCutoff):
                print '%s-%s\t%f\t%f\t%f' % (lineObj.chrom(), lineObj.loc(), lineObj.ao(), lineObj.dp(), lineObj.af())

if __name__ == '__main__':
    inFile, afCutoff, aoCutoff = runArgParse()
    readFile(inFile, afCutoff, aoCutoff)



