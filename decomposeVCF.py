# this takes care of the problem in freebayes where variants are output
# not left aligned and not parsimonious
# for ref: http://genome.sph.umich.edu/wiki/Variant_Normalization
from os import system

def decompose(vcfOut):
    decomposeOut = vcfOut.strip('.vcf') + 'Decomposed.vcf'
    blockDecomposedOut = vcfOut.strip('.vcf') + 'BlockDecomposed.vcf'
    system('vt decompose %s > %s' % (vcfOut, decomposeOut))
    system('vt decompose_blocksub %s > %s' % (decomposeOut, blockDecomposedOut))

    return blockDecomposedOut
