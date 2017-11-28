#!/usr/bin/env python

# this takes a line from a vcf file and outputs sequence of specified length
# flanking mutated locus
def getSequence(vcfLine, flankLength):
    from parseline import VCFObj
    from subprocess import check_output, STDOUT

    vcfObj = VCFObj(vcfLine)

    low = int(vcfObj.location) - flankLength
    high = int(vcfObj.location) + flankLength
    temp = check_output('wget -qO- http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=%s:%s,%s' % (vcfObj.chrom,low,high), stderr=STDOUT, shell=True)

    finalSeq = ''
    for line in temp.split('\n'):
        if '<' not in line:
            finalSeq += line

    return finalSeq

# this provides similar functionality to getSequence, but instead of using the
# USCS genome browser for sequence, it uses a reference genome
def getRefSequence(vcfLine, flankLength, ref):
    from parseline import VCFObj
    from subprocess import check_output, STDOUT
    from string import upper

    vcfObj = VCFObj(vcfLine)

    low = int(vcfObj.location) - flankLength
    high = int(vcfObj.location) + flankLength
    temp = check_output('samtools faidx %s %s:%s-%s' % (ref, vcfObj.chrom, low, high), stderr=STDOUT, shell=True)

    finalSeq = ''
    for line in temp.split('\n'):
        if '>' not in line:
            finalSeq += line

    finalSeq = finalSeq.upper()
    return finalSeq




