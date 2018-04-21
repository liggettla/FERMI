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
parser.add_argument('--justoncosites', '-j', action='store_true', help='This will just plot oncogenic variants.')
parser.add_argument('--pythonplotting', '-n', action='store_true', help='This will trigger use of the python plotting algorithm.')


args = parser.parse_args()

# specifiy which plotting algorithm to use
if args.pythonplotting:
    pythonplotting = True
else:
    pythonplotting = False

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

if args.justoncosites:
    justoncosites = True
else:
    justoncosites = False

#################
# Specify Files #
#################
inputDir = args.indir
outputDir = args.outdir
#principle = inputDir + '/' + args.principle + '/' +'total_filtered.vcf'
#principle = inputDir + '/' + args.principle + '/' +'onlyProbedRegions.vcf'
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

    location = '%s-%s-%s-%s' % (chrom, str(loc), str(WT), str(var))
    loc = int(loc)
    return location, AFNum, WT, var, loc

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
def buildAverageStructure(samples, regions, oncosites, justoncosites):
    from Bio.Seq import Seq
    tempData = {}
    for i in samples:
        target = open(inputDir + '/' + i + '/' +'onlyProbedRegions.vcf', 'r')
        for line in target:
            if '#' not in line and 'chr' in line: # skip the info
                # Ex: loc = chr1-1234-C-A
                loc, AFNum, WT, var, location = parseLine(line)
                if args.multiplier:
                    AFNum = AFNum * args.multiplier

                # decide if in acceptable probe regions; coding vs template
                withinProbes = False
                if not justoncosites:
                    for i in regions:
                        start=int(regions[i][0])
                        end=int(regions[i][1])
                        if location > start and location < end:
                            withinProbes = True

                    if withinProbes:
                        # should germline/variant types be included?
                        # and either combine complementary bases or not
                        if args.combinecomplements:
                            matchingVariants = WT in germline and var in variant or str(Seq(WT).complement()) in germline and str(Seq(var).complement()) in variant
                        else:
                            matchingVariants = WT in germline and var in variant

                        if matchingVariants:
                            # decide if variant is unique or not
                            if rareEnough(AFNum):
                                if loc in tempData:
                                    tempData[loc]['vaf'].append(AFNum)
                                else:
                                    tempData[loc] = {'vaf':[AFNum]}
                else:
                    if str(location) in oncosites:

                        if args.combinecomplements:
                            matchingVariants = WT in germline and var in variant or str(Seq(WT).complement()) in germline and str(Seq(var).complement()) in variant
                        else:
                            matchingVariants = WT in germline and var in variant

                        if matchingVariants:
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
def buildPrincipleStructure(principle, regions, oncosites, justoncosites):
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
    from os import system
    #import pdb; pdb.set_trace()

    # plot the results
    command = 'Rscript plotvafRepeatability.R'
    if displayPlot:
        system(command)

    # move the results
    command = 'mv outputFile %s' % (outputFile)
    system(command)
    #command = 'mv output1.jpg %s' % (plotFile1)
    #system(command)
    command = 'mv output2.jpg %s' % (plotFile2)
    if displayPlot:
        system(command)

    # display the plots
    #command = 'eog vafRepeatabilityRegression.jpg'
    #system(command)
    command = 'eog vafRepeatabilityNoRegression.jpg'
    if displayPlot:
        system(command)

##################
# Run the Script #
##################
if __name__ == '__main__':
    oncosites = ['5073770','7577539','7577119','115256529','115258747','115258744','534287','534288','534289','25398284','25380275','106197266','106197267','106197268','106197269','106155172','106155173','106155174','25457242','25457243','209113112','209113113','90631934','90631838','48649700']
    if args.onlyCoding:
        strand = 'Coding'
        regions = {'TIIIb':['223190673','223190820'],'TET2-1':['106197237','106197405'],'TET2-2':['106155137','106155275'],'TIIId':['110541172','110541302'],'JAK2':['5073733','5073887'],'TIIIl':['2593889','2594074'],'TIIIn':['92527052','92527176'],'TIIIq':['85949137','85949299'],'GATA1':['48649667','48649849']}

    elif args.onlyTemplate:
        strand = 'Template'
        regions = {'TIIIa':['115227813','115227978'],'NRAS-1':['115256496','115256680'],'NRAS-2':['115258713','115258897'],'DNMT3a':['25457211','25457364'],'IDH1':['209113077','209113239'],'SF3B1':['198266803','198266967'],'TIIIc':['229041101','229041289'],'TIIIk':['2389983','2390171'],'TIIIl':['2593889','2594074'],'TIIIm':['11486596','11486728'],'HRAS':['534258','534385'],'KRAS-1':['25398247','25398415'],'KRAS-2':['25380242','25380368'],'IDH2':['90631809','90631969'],'p53-1':['7577504','7577635'],'p53-2':['7578369','7578544'],'p53-3':['7577084','7577214']}

    else:
        regions = {'TIIIa':['115227814','115227978'],'NRAS-1':['115256496','115256680'],'NRAS-2':['115258713','115258897'],'DNMT3a':['25457211','25457364'],'IDH1':['209113077','209113239'],'SF3B1':['198266803','198266967'],'TIIIb':['223190674','223190820'],'TIIIc':['229041101','229041289'],'TET2-1':['106197237','106197405'],'TET2-2':['106155137','106155275'],'TIIId':['110541172','110541302'],'TIIIe':['112997214','112997386'],'TIIIf':['121167756','121167884'],'TIIIg':['123547743','123547901'],'TIIIh':['124428637','124428767'],'JAK2':['5073733','5073887'],'TIIIj':['2126256','2126420'],'TIIIk':['2389983','2390171'],'TIIIl':['2593889','2594074'],'TIIIm':['11486596','11486728'],'HRAS':['534258','534385'],'KRAS-1':['25398247','25398415'],'KRAS-2':['25380242','25380368'],'TIIIn':['92527052','92527176'],'IDH2':['90631809','90631969'],'TIIIo':['73379656','73379832'],'TIIIp':['82455026','82455164'],'TIIIq':['85949137','85949299'],'p53-1':['7577504','7577635'],'p53-2':['7578369','7578544'],'p53-3':['7577084','7577214'],'GATA1':['48649667','48649849']}

    avgData = buildAverageStructure(samples, regions, oncosites, justoncosites)

    # principle can either be a single sample or multiple samples
    # to save time this will be conditionally handled by existing methods
    if len(args.principle) == 1:
        principle = inputDir + '/' + args.principle[0] + '/' +'onlyProbedRegions.vcf'
        principleData = buildPrincipleStructure(principle, regions, oncosites, justoncosites)

    else:
        principleData = buildAverageStructure(args.principle, regions, oncosites, justoncosites)

    # this will use the original R plotting algorithm
    if not pythonplotting:
        outputData(commonVars, avgData, principleData)

        # output plot if requested
        if args.displayplot:
            plotAndDisplay(outputFile, plotFile1, plotFile2, True)
        else:
            plotAndDisplay(outputFile, plotFile1, plotFile2, False)

        # this now also reports pearsons correlation coefficient
        from revisedComputeRSquared import getRSquared
        r2, p = getRSquared()
        print('R-Squared = %f' % (r2))
        print('Pearson Coefficient = %f' % (p[0]))
        print('p-value = %f' % (p[1]))

        if args.plotonchrom:
            from os import system
            system('python plotOnChromosome.py')

    # this will use the newer and hopefully superior python plotting algorithm
    if pythonplotting:
        from  plotvafRepeatability import builddf
        builddf(commonVars, avgData, principleData)


