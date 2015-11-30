#The purpose of this script is to make use of both UMIs found in paired-end
#reads. The problem is that the HiSeq chemistry only runs for 100 cycles
#which means that unless stitching works well, both UMIs are never found
#in a single read file. This script overcomes this problem by reading the
#UMI sequences from R2 and concatentating them onto the corresponding reads
#in R1. These sequences can then be binned based on their UMI sequences
#while ignoring the problems that stitching can introduce.

import itertools
from Bio.Seq import Seq

#use entire read from read1, and the last 6bp of read2; concatenate together
#and write to output
def concatenateUMI(read1, read2, output):
    position = 1
    plus = '+'
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
                    r2UMI = revComp[-6:]

                    position += 1
                elif position == 3:
                    position += 1
                elif position == 4:
                    quality = line1.rstrip('\n')
                    r2UMIqual = line2[0:6]
                    position = 1
                    target.write(header + '\n' + seq1 + r2UMI + '\n' + plus + '\n' + quality + r2UMIqual + '\n')
