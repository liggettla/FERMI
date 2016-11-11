#!/usr/bin/python
# the purpose of this script is to analyse blood samples that were prepped
# twice to understand if variants that are unique to the samples when compared
# to an average are still unique when compared to each other.

############
# Argparse #
############
def runArgParse():
    import argparse
    from numpy import mean

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', '-i', type=str, required=True, help='Input directory containing the vcf files to be analyzed: /dir')
    parser.add_argument('--outdir', '-o', type=str, required=True, help='Output directory for plots: /dir')
    parser.add_argument('--averageSamples', '-c', type=str, nargs='*', required=True, help='List of samples to be averaged and compared to the principle samples. Ex: A1-R1')
    parser.add_argument('--sampleOne', '-a', type=str, nargs='*', required=True, help='The first of the two samples to be compared to average and to its technical replicate.')
    parser.add_argument('--sampleTwo', '-b', type=str, nargs='*', required=True, help='The second of the two samples to be compared to average and to its technical replicate.')

    args = parser.parse_args()

    inputDir = args.indir
    outputDir = args.outdir
    avgSamples = args.averageSamples # this is a list
    sampleOne = args.sampleOne
    sampleTwo = args.sampleTwo

    return inputDir, outputDir, avgSamples, sampleOne, sampleTwo

if __name__ == '__main__':
    from buildDF import buildDF

    inputDir, outputDir, avgSamples, sampleOne, sampleTwo = runArgParse()

    dfAvg = buildDF(inputDir, avgSamples)
    df1 = buildDF(inputDir, sampleOne)
    df2 = buildDF(inputDir, sampleTwo)
