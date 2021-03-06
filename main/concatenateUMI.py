#####################
#Concatate R1 and R2#
#####################
# The purpose of this script is to make use of both UMIs found in paired-end
# reads. The problem is that the HiSeq chemistry only runs for 100 cycles
# which means that unless stitching works well, both UMIs are never found
# in a single read file. This script overcomes this problem by reading the
# UMI sequences from R2 and concatentating them onto the corresponding reads
# in R1. These sequences can then be binned based on their UMI sequences
# while ignoring the problems that stitching can introduce.
# attach 3' UMI from R2 onto R1 read
# this is a necessary step to process reads with 100-150 cycle chemistry
# 200 cycle chemistry should make this unnecessary

import itertools
from Bio.Seq import Seq

#use entire read from read1, and the last 6bp of read2; concatenate together
#and write to output
def concatenateUMI(read1, read2, output):
    position = 1
    plus = '+'
    #write results to a new fastq file
    target = open(output, 'w')

    with open(read1, 'r') as R1:
        with open(read2, 'r') as R2:

            #in python2 zip reads everything into memory
            #itertools.izip is iterative, and reads line by line
            for line1, line2 in itertools.izip(R1, R2):
                if position == 1:
                    header = line1.rstrip('\n')
                    position += 1
                elif position == 2:
                    seq1 = line1.rstrip('\n')

                    #It is necessary to first strip, or \n will count as a char
                    seq2 = line2.rstrip('\n')
                    seq2Obj = Seq(seq2)

                    #use read1 strand
                    revComp = str(seq2Obj.reverse_complement())
                    r2UMI = revComp[-11:]

                    position += 1
                elif position == 3:
                    position += 1
                #retain quality scores for R2 UMI
                elif position == 4:
                    quality = line1.rstrip('\n')
                    r2UMIqual = line2[0:6]
                    position = 1
                    target.write(header + '\n' + seq1 + r2UMI + '\n' + plus + '\n' + quality + r2UMIqual + '\n')

def duplexConcatenate(read1, read2, output, realvsmock, overlap=80):
    position = 1
    plus = '+'

    #write results to a new fastq file
    target = open(output, 'w')

    with open(read1, 'r') as R1:
        with open(read2, 'r') as R2:

            #in python2 zip reads everything into memory
            #itertools.izip is iterative, and reads line by line
            for line1, line2 in itertools.izip(R1, R2):
                if position == 1:
                    header = line1.rstrip('\n')
                    position += 1
                elif position == 2:
                    seq1 = line1.rstrip('\n')
                    seq2 = str(Seq(line2.rstrip('\n')).reverse_complement())
                    r1UMI = seq1[:6]
                    r2UMI = seq2.rstrip('\n')[-6:]

                    r1 = seq1[6:]
                    r2 = seq2[:-6]
                    r1,r2 = trimReads(r1,r2)

                    finalRead = duplexCollapse(r1, r2, r1UMI, r2UMI, realvsmock)

                    position += 1
                elif position == 3:
                    position += 1
                elif position == 4:
                    quality = line1.rstrip('\n')[:len(finalRead)]
                    position = 1

                    # only write reads that are of significant length
                    # this should only be necessary when duplex collapsing
                    if len(r1) > 5 and len(r2) > 5:
                        target.write(header + '\n' + finalRead + '\n' + plus + '\n' + quality + '\n')

def duplexCollapse(r1, r2, r1UMI, r2UMI, realvsmock):
    import numpy as np
    finalRead = ''
    finalRead += r1UMI
    # this is crude, but if r1 and r2 are not of the same length just use r1
    # it would probably be better to throw it out but for now this way is easier
    if len(r1) == len(r2):
        for i in np.arange(len(r1)):
            if r1[i] == r2[i]:
                finalRead += r1[i]
            elif r1[i] != r2[i]:
                # this will either fix the error or use the base from r1 in order to understand the
                # effect duplex collapsing is having
                if realvsmock:
                    finalRead += 'N'
                else:
                    finalRead += r1[i]

        finalRead += r2UMI
        return finalRead

    else:
        return r1

def trimReads(r1,r2):
    from Levenshtein import distance

    # some of the captures are less than 150bp causing r2 to extend 5' of r1
    # so cleave off the most that could ever overhang from my captures and now
    # r1 should always align more 5' than r2
    r1 = r1[:-20]
    r2 = r2[20:]

    r1_5prime = r1[:10] # first 10bp
    r1_3prime = r1[-10:]# last 10bp
    # reversed for r2
    r2_5prime = r2[-10:]
    r2_3prime = r2[:10]

    while distance(r1_5prime, r2_3prime) > 0 and len(r1)>5: # allow one error and trim until read almost gone
        r1 = r1[1:] # cut off first base
        r1_5prime = r1[:10] # first 10bp

    while distance(r2_5prime, r1_3prime) > 0 and len(r2)>5: # allow one error
        r2 = r2[:-1] # cut off last base
        r2_5prime = r2[-10:] # first 10bp

    return r1, r2

















