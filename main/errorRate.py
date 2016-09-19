'''
This script estimates the overall base reporting error rate
including pcr errors, sequecing errors, and other error sources.
This estimationg is done by totaling the number of bases that are not
used as the consensus base and dividing that by consensus base observation
times. This is done for every analyzed base locus to derive an average
error rate.

This estimation may be slightly high as no method for dealing with two
different captures and the same UMI has been explicitly written yet.
These events will likely trend to 50% error rates
'''
from numpy import mean

def calcErrorRate(a,t,g,c):
    bases = [a,t,g,c]
    consensus = float(max(bases))
    totalReads = float(sum(bases))

    bases.remove(max(bases))
    errors = float(sum(bases))

    errorRate = errors / totalReads * 100
    return errorRate
