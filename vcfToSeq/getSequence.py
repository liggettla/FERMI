#!/usr/bin/env python

# this takes a line from a vcf file and outputs sequence of specified length
# flanking mutated locus
def getSequence(vcfLine, flankLength):
    from parseline import parseLine
    from subprocess import check_output, STDOUT

    location, AFNum, WT, var, loc, chrom = parseLine(vcfLine)

    low = int(loc) - flankLength
    high = int(loc) + flankLength
    temp = check_output('wget -qO- http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=%s:%s,%s' % (chrom,low,high), stderr=STDOUT, shell=True)

    finalSeq = ''
    for line in temp.split('\n'):
        if '<' not in line:
            finalSeq += line

    return finalSeq
