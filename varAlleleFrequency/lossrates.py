#!/usr/bin/env python

############
# Argparse #
############
def runArgparse():
    import argparse
    from numpy import mean

    parser = argparse.ArgumentParser()
    parser.add_argument('--inFile', '-i', type=str, help='Specifies the input vcf file.')

    args = parser.parse_args()
    inFile = args.inFile

    return inFile

def parseFile(inFile):
    from parseline import parseLine
    target = open(inFile, 'r')

    baseCounts = {'purine':0, 'pyrimidine':0}
    purines = ['A', 'G']
    pyrimidines = ['C', 'T']
    for line in target:
        location, AFNum, WT, var, loc = parseLine(line)
        baseType = 'purine' if var in purines else 'pyrimidine'
        baseCounts[baseType] += 1

    print baseCounts
    return baseCounts

def displayCounts(baseCounts):
    for i in baseCounts:
        print i, baseCounts[i]

if __name__ == '__main__':
    inFile = runArgparse()
    baseCounts = parseFile(inFile)
    displayCounts(baseCounts)
