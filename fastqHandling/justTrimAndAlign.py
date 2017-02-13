#!/usr/bin/env python
# The purpose of this script is to just trim and align reads
# in order to compare how much of a difference the processing
# is making.

def argParse():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', '-i', type=str, help='Location of input file.')
    parser.add_argument('--filename', '-f', type=str, help='The name of the file to be processed.')

    args = parser.parse_args()
    inDir = args.indir
    filename = args.filename
    return inDir, filename

def trimUMI(inDir, inFile):
    target = open(inFile, 'r')
    outTarget = open(inDir + '/noUMI', 'w')
    position = 1

    for line in target:
        if position == 1:
            header = line.rstrip('\n')
            position += 1
        elif position == 2:
            #Assumes UMI is flanking first and last 6bp of read
            umi_seq = line[0:11]+line[-6:] #Abs dist from start/end compatible with miSeq/hiSeq
            umi_seq = umi_seq.rstrip('\n')
            read_seq = line[6:25]
            position += 1
        elif position == 3:
            position += 1
        elif position == 4:
            quality = line.rstrip('\n')
            position = 1

    readSeq = line[6:-25]

if __name__ == '__main__':
    inDir, filename = argParse()
    inFile = inDir + '/' + filename
    trimUMI(inDir, inFile)
