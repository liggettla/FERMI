#!/usr/bin/python
# the purpose of this script is to compare two repeatedly sequenced samples originating
# from the same DNA sample to understand how repeatable the VAF is for a given mutation
# in each of the analysis files
# for example in 1.vcf VAF of a variant might be 0.5 but in 2.vcf it is 0.2

############
# Argparse #
############
import argparse
from numpy import mean

parser = argparse.ArgumentParser()
parser.add_argument('--indir', '-i', type=str, required=True, help='Input directory containing the vcf files to be analyzed: /dir')
parser.add_argument('--outdir', '-o', type=str, required=True, help='Output directory for plots: /dir')
parser.add_argument('--principle', '-p', type=str, required=True, help='This is the principle sample being compared to an averaged set of other samples. Ex: A1-R1')
parser.add_argument('--samples', '-s', type=str, nargs='*', required=True, help='List of samples to be averaged and compared to the principle sample. Ex: A1-R1')
parser.add_argument('--rarevars', '-r', type=float, help='This can be set to cutoff the data at a certain allele frequency and only include variants below a particular frequency like 0.03 or 0.003.')
parser.add_argument('--commonVars', '-c', action='store_true', help='This will only plot variants that are found in both samples and ignore those variants that are only found in one of the samples.')

args = parser.parse_args()

inputDir = args.indir
outputDir = args.outdir
#principle = inputDir + '/' + args.principle + '/' +'total_filtered.vcf'
principle = inputDir + '/' + args.principle + '/' +'onlyProbedRegions.vcf'
samples = args.samples # this is a list
commonVars = args.commonVars
cutoff = args.rarevars

# define output files
outputFile = outputDir + '/vafRepeatability.txt'
plotFile1 = outputDir + '/vafRepeatabilityRegression.jpg'
plotFile2 = outputDir + '/vafRepeatabilityNoRegression.jpg'

################
# Check Rarity #
################
# check if rare variants are expected
# and if so, if variant is rare enough
def rareEnough(AFNum):
    if not cutoff:
        return True
    if cutoff and AFNum <= cutoff:
        return True
    else:
        return False

##############
# Parse Line #
##############
def parseLine(i):
    chrom = str(i.split('\t')[0])
    loc = str(i.split('\t')[1])
    AO = i.split(';')[5]
    DP = i.split(';')[7]
    var = i.split('\t')[4] # the observed bp change
    WT = i.split('\t')[3]

    AONum = float(AO.split(',')[0][3:])
    DPNum = float(DP.split(',')[0][3:])
    AFNum = AONum / DPNum

    location = '%s-%s-%s' % (chrom, str(loc), str(var))
    return location, AFNum

###################
# Find Avg AFNums #
###################
# computes the average AFNum for each unique variant
# and returns averaged data structure
def takeAverage(tempData):
    avgData = {}
    for loc in tempData:
        x = mean(tempData[loc]['vaf'])
        avgData[loc] = x

    return avgData

############################
# Build Avg Data Structure #
############################
def buildAverageStructure(samples):
    tempData = {}
    for i in samples:
        target = open(inputDir + '/' + i + '/' +'total_filtered.vcf', 'r')
        for line in target:
            if '#' not in line and 'chr' in line: # skip the info
                # Ex: loc = chr1:1234:A
                loc, AFNum = parseLine(line)

                # decide if variant is unique or not
                if rareEnough(AFNum):
                    if loc in tempData:
                        tempData[loc]['vaf'].append(AFNum)
                    else:
                        tempData[loc] = {'vaf':[AFNum]}

    # average all of the AFNum values
    avgData = takeAverage(tempData)
    return avgData

#############################
# Build Principle Structure #
#############################
# Builds dictionary of the principle data
# to be compared to the average data
def buildPrincipleStructure(principle):
    principleData = {}
    target = open(principle, 'r')
    for i in target:
        if '#' not in i and 'chr' in i: # skip the info
            loc, AFNum = parseLine(i)
            # decide if variant is unique or not
            if rareEnough(AFNum):
                principleData[loc] = AFNum

    return principleData

#########################
# Output Plottable Data #
#########################
# writes data to temporary file that is then read in by R for plotting
def outputData(commonVars, avgData, principleData):
    output = open('outputFile', 'w')
    output.write('Sample1\tSample2\tIdentity\n')

    # only output variants observed in both samples
    if commonVars:
        for i in avgData:
            if i in principleData:
                output.write('%s\t%s\t%s\n' % (principleData[i], avgData[i], i))

    # output all variants including those not observed in both samples
    else:
        for i in avgData:
            # write overlapping variants
            if i in principleData:
                output.write('%s\t%s\t%s\n' % (principleData[i], avgData[i], i))
            # write variants found only in avgData
            else:
                output.write('%s\t%s\t%s\n' % (0, avgData[i], i))

        # write variants found only in principleData
        for i in principleData:
            if i not in avgData:
                output.write('%s\t%s\t%s\n' % (principleData[i], 0, i))

    output.close()

####################
# Plot and Display #
####################
def plotAndDisplay(outputFile, plotFile1, plotFile2):
    print outputFile, plotFile1, plotFile2
    from os import system
    #import pdb; pdb.set_trace()

    # plot the results
    command = 'Rscript plotvafRepeatability.R'
    system(command)

    # move the results
    command = 'mv outputFile %s' % (outputFile)
    system(command)
    command = 'mv output1.jpg %s' % (plotFile1)
    system(command)
    command = 'mv output2.jpg %s' % (plotFile2)
    system(command)

    # display the plots
    #command = 'eog vafRepeatabilityRegression.jpg'
    #system(command)
    command = 'eog vafRepeatabilityNoRegression.jpg'
    system(command)

##################
# Run the Script #
##################
avgData = buildAverageStructure(samples)
principleData = buildPrincipleStructure(principle)
outputData(commonVars, avgData, principleData)
plotAndDisplay(outputFile, plotFile1, plotFile2)
