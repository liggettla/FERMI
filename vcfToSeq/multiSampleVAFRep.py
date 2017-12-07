#!/usr/bin/python
# the purpose of this script is to compare two repeatedly sequenced samples originating
# from the same DNA sample to understand how repeatable the VAF is for a given mutation
# in each of the analysis files
# for example in 1.vcf VAF of a variant might be 0.5 but in 2.vcf it is 0.2

############
# Argparse #
############
def runArgparse():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', '-i', type=str, required=True, help='Input directory containing the vcf files to be analyzed: /dir')
    parser.add_argument('--outdir', '-o', type=str, required=True, help='Output directory for plots: /dir')
    parser.add_argument('--principle', '-p', type=str, nargs='*', required=True, help='These are the principle samples being compared to an averaged set of other samples. Ex: A1-R1')
    parser.add_argument('--samples', '-s', type=str, nargs='*', required=True, help='List of samples to be averaged and compared to the principle sample. Ex: A1-R1')
    parser.add_argument('--rarevars', '-r', type=float, help='This can be set to cutoff the data at a certain allele frequency and only include variants below a particular frequency like 0.03 or 0.003.')
    parser.add_argument('--commonVars', '-c', action='store_true', help='This will only plot variants that are found in both samples and ignore those variants that are only found in one of the samples.')
    parser.add_argument('--germline', '-g', type=str, nargs='*', help='Only output those variants that changed from these bases.')
    parser.add_argument('--variant', '-v', type=str, nargs='*', help='Only output those variants that change to these bases.')
    parser.add_argument('--displayplot', '-d', action='store_true', help='This will trigger the displaying of VAF plot.')
    parser.add_argument('--multiplier', '-m', type=float, help='This specifies a multiplier to artificially increase all the VAFs in the principle sample by a set multiplier allowing for comparison of shifting populations.')
    parser.add_argument('--plotonchrom', '-z', action='store_true', help='This will output vafs of the data shown on the y-axis of the vaf comparison plot (typically the average samples) along chromosomal distances to understand hot and cold regions of the chromosome.')
    parser.add_argument('--combinecomplements', '-a', action='store_true', help='This will combine the complement of base pairs into a single plot, ie if C-T variants are asked for, both C-T and G-A variants will be output.')
    parser.add_argument('--onlyCoding', '-x', action='store_true', help='This will only use variants coming from coding strand probes.')
    parser.add_argument('--onlyTemplate', '-y', action='store_true', help='This will only use variants coming from template strand probes.')
    parser.add_argument('--tripletcontext', '-t', type=str, help='If only one triplet context should be analyzed, this can be used to specify that context.')
    parser.add_argument('--referencegenome', '-f', type=str, help='Point to the reference genome to be used if flanking triplet context is to be used in analysis')
    parser.add_argument('--loadpreviousdata', '-l', action='store_true', help='Analysis takes a long time if triplet context is being considered (up to 1min per vcf file), so this can be used to load previously processed data by pointing to previous output pickle file.')

    args = parser.parse_args()

    if args.displayplot:
        displayplot = True
    else:
        displayplot = False

    if args.plotonchrom:
        plotonchrom = True
    else:
        plotonchrom = False

    if args.multiplier:
        multiplier = args.multiplier
    else:
        multiplier = False

    if args.combinecomplements:
        combinecomplements = True
    else:
        combinecomplements = False

    principlefile = args.principle

    #####################
    # Specify Mutations #
    #####################
    if args.germline:
        germline = args.germline
    else:
        germline = ['A','T','G','C']

    if args.variant:
        variant = args.variant
    else:
        variant = ['A','T','G','C']

    if args.tripletcontext:
        from string import upper
        triple = upper(args.tripletcontext)
    else:
        triple = ''

    if args.referencegenome:
        ref = args.referencegenome
    else:
        ref = ''

    if args.onlyCoding:
        strand = 'Coding'
    elif args.onlyTemplate:
        strand = 'Template'
    else:
        strand = 'All'

    #################
    # Specify Files #
    #################
    inputDir = args.indir
    outputDir = args.outdir
    samples = args.samples # this is a list
    commonVars = args.commonVars
    cutoff = args.rarevars

    # define output files
    outputFile = outputDir + '/vafRepeatability.txt'
    plotFile1 = outputDir + '/vafRepeatabilityRegression.jpg'
    plotFile2 = outputDir + '/vafRepeatabilityNoRegression.jpg'

    # use previous data or not
    if args.loadpreviousdata:
        previousData = True
    else:
        previousData = False

    return inputDir, outputDir, samples, commonVars, cutoff, outputFile, plotFile1, plotFile2, germline, variant, triple, ref, strand, previousData, principlefile, displayplot, plotonchrom, multiplier, combinecomplements

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

    location = '%s-%s-%s-%s' % (chrom, str(loc), str(WT), str(var))
    loc = int(loc)
    return location, AFNum, WT, var, loc

###################
# Find Avg AFNums #
###################
# computes the average AFNum for each unique variant
# and returns averaged data structure
def takeAverage(tempData):
    from numpy import mean

    # tempData = {(1500, 'CAC'):[0.75,0.35,0.68]}
    avgData = {}
    for loc in tempData:
        avgData[loc] = mean(tempData[loc])

    return avgData

############################
# Build Avg Data Structure #
############################
def buildAverageStructure(samples, regions, multiplier, combinecomplements):
    from Bio.Seq import Seq
    from getSequence import getRefSequence
    from collections import defaultdict

    tempData = {}
    holdingData = defaultdict(list)
    for i in samples:
        target = open(inputDir + '/' + i + '/' +'onlyProbedRegions.vcf', 'r')
        for line in target:
            if '#' not in line and 'chr' in line: # skip the info
                # Ex: loc = chr1-1234-C-A
                loc, AFNum, WT, var, location = parseLine(line)
                if multiplier:
                    AFNum = AFNum * multiplier

                '''
                # decide if in acceptable probe regions; coding vs template
                withinProbes = False
                for i in regions:
                    start=int(regions[i][0])
                    end=int(regions[i][1])
                    if location > start and location < end:
                        withinProbes = True

                if withinProbes:
                    # should germline/variant types be included?
                    # and either combine complementary bases or not
                    if combinecomplements:
                        matchingVariants = WT in germline and var in variant or str(Seq(WT).complement()) in germline and str(Seq(var).complement()) in variant
                    else:
                        matchingVariants = WT in germline and var in variant
                '''

                # eliminate SNPs
                if AFNum < 0.1:
                    # get flanking sequence
                    seq = getRefSequence(line, 1, ref)

                    # tempData = {(1500, 'CAC'):[0.75,0.35,0.68]}
                    holdingData[(loc, seq)].append(AFNum)

    # average all of the AFNum values
    avgData = takeAverage(holdingData)
    return avgData

#############################
# Build Principle Structure #
#############################
# Builds dictionary of the principle data
# to be compared to the average data
def buildPrincipleStructure(principle, regions):
    from Bio.Seq import Seq
    principleData = {}
    target = open(principle, 'r')
    for i in target:
        if '#' not in i and 'chr' in i: # skip the info
            loc, AFNum, WT, var, location = parseLine(i)

            # should germline/variant types be included?
            if WT in germline and var in variant or str(Seq(WT).complement()) in germline and str(Seq(var).complement()) in variant:

                # decide if in acceptable probe regions; coding vs template
                withinProbes = False
                for i in regions:
                    start=int(regions[i][0])
                    end=int(regions[i][1])
                    if location > start and location < end:
                        withinProbes = True
                if withinProbes:
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
def plotAndDisplay(outputFile, plotFile1, plotFile2, displayPlot):
    from subprocess import call

    # plot the results
    command = 'Rscript plotvafRepeatability.R'
    if displayPlot:
        call(command, shell=True)

    # move the results
    command = 'mv outputFile %s' % (outputFile)
    call(command, shell=True)
    command = 'mv output2.jpg %s' % (plotFile2)
    if displayPlot:
        call(command, shell=True)

    # display the plots
    command = 'eog vafRepeatabilityNoRegression.jpg'
    if displayPlot:
        call(command, shell=True)

def getPreviousData():
    import pickle
    p = open('principle.pkl', 'rb')
    s = open('samples.pkl', 'rb')
    principleData = pickle.load(p)
    avgData = pickle.load(s)
    p.close()
    s.close()

    return principleData, avgData

def getNewData(samples, regions, inputDir, principlefile, multiplier, combinecomplements):
    # y-axis
    avgData = buildAverageStructure(samples, regions, multiplier, combinecomplements)

    # x-axis
    if len(principlefile) == 1:
        principle = inputDir + '/' + principlefile[0] + '/' +'onlyProbedRegions.vcf'
        principleData = buildPrincipleStructure(principle, regions)

    else:
        principleData = buildAverageStructure(principlefile, regions, multiplier, combinecomplements)

    return avgData, principleData

def writeNewData(avgData, principleData):
    import pickle
    p = open('principle.pkl', 'wb')
    s = open('samples.pkl', 'wb')
    pickle.dump(principleData, p)
    pickle.dump(avgData, s)
    p.close()
    s.close()

def calculateStats():
    from revisedComputeRSquared import getRSquared
    r2, p = getRSquared()
    print('R-Squared = %f' % (r2))
    print('Pearson Coefficient = %f' % (p[0]))
    print('p-value = %f' % (p[1]))

# this will eliminate any triplets that were not requested
def filterTriplets(avgData, principleData, triple, germline, variant, combinecomplements):
    from numpy import mean
    from Bio.Seq import Seq
    a = {}
    p = {}

    g = []
    v = []
    if combinecomplements:
        for i in germline:
            g.append(i)
            g.append(str(Seq(i).complement()))
        for i in variant:
            v.append(i)
            v.append(str(Seq(i).complement()))


    for i in avgData:
        if i[1] == triple or str(Seq(i[1]).reverse_complement()) == triple:
            if i[0].split('-')[2] in g and i[0].split('-')[3] in v:
                a[i[0]] = mean(avgData[i])

    for i in principleData:
        if i[1] == triple or str(Seq(i[1]).reverse_complement()) == triple:
            if i[0].split('-')[2] in g and i[0].split('-')[3] in v:
                p[i[0]] = mean(principleData[i])

    return a, p

##################
# Run the Script #
##################
if __name__ == '__main__':
    print('Parsing Args...')
    inputDir, outputDir, samples, commonVars, cutoff, outputFile, plotFile1, plotFile2, germline, variant, triple, ref, strand, previousData, principlefile, displayplot, plotonchrom, multiplier, combinecomplements = runArgparse()

    # specify which probes to use ie: coding or template
    print('Specify probes...')
    from personalRegions import probeLocation
    regions = probeLocation(strand)

    # use previous data structures
    if previousData:
        print('Loading Previous Data...')
        principleData, avgData = getPreviousData()

    # build new data structures
    else:
        print('Building New Data Structures (Takes up to 1min per input file)...')
        avgData, principleData = getNewData(samples, regions, inputDir, principlefile, multiplier, combinecomplements)

        print('Writing New Data...')
        writeNewData(avgData, principleData)

    # filter requested triplets
    print('Filtering Triplets...')
    if len(triple) > 1:
        avgData, principleData = filterTriplets(avgData, principleData, triple, germline, variant, combinecomplements)

    # write temp file for R plotting
    print('Feed Data Into R...')
    outputData(commonVars, avgData, principleData)

    # output plot if requested
    if displayplot:
        print('Plotting VAF Comparison...')
        plotAndDisplay(outputFile, plotFile1, plotFile2, True)
    else:
        print('Calculating Comparison...')
        plotAndDisplay(outputFile, plotFile1, plotFile2, False)

    # this now also reports pearsons correlation coefficient
    print('Calculting Statistics...')
    calculateStats()

    # display confetti plot if requested
    if plotonchrom:
        from os import system
        system('python plotOnChromosome.py')




