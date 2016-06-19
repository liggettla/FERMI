###############
#Call Variants#
###############
#calls variants only if they are at least 0.00001% of the calls 1/10^6
#this should thus output every variant found in the alignment
#by default freebayes uses 0.2 (20%)

from os import system

def callVar(REF, bamOut, vcfOut):

    system("%s -X -F 0.0000001 --fasta-reference %s %s > %s" % (freebayes, REF, bamOut, vcfOut))

#Using pooled continuous makes frequency based calls without using number of samples as input
#This might help in adjusting for different copy numbers without knowing the exact number of individual captures
#in a given reaction
#freebayes -f %s -F 0.0000001 -C 1 --pooled-continuous %s > %s % (REF, bamOut, vcfOut))

