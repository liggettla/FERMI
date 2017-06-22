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

def duplexConcatenate(read1, read2, output, overlap=80):
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
                    r2UMI = seq2[-6:]
                    r1 = seq1[6:]
                    r2 = seq2[:-6]
                    r1 = r1[-overlap:]
                    r2 = r2[:overlap]
                    print r1, r2, r1UMI, r2UMI

                    finalRead = duplexCollapse(r1, r2, r1UMI, r2UMI)

                    position += 1
                elif position == 3:
                    position += 1
                #retain quality scores for R2 UMI
                elif position == 4:
                    quality = line1.rstrip('\n')[:len(finalRead)]
                    position = 1
                    target.write(header + '\n' + finalRead + '\n' + plus + '\n' + quality + '\n')

def duplexCollapse(r1, r2, r1UMI, r2UMI):
    import numpy as np
    finalRead = ''
    finalRead += r1UMI
    for i in np.arange(len(r1)):
        if r1[i] == r2[i]:
            finalRead += r1[i]
        elif r1[i] != r2[i]:
            finalRead += 'N'

    finalRead += r2UMI
    return finalRead






