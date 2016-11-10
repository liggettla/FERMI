#!/usr/bin/python
# the purpose of this script is to build different dataframes that parse
# through and associate different details from vcf files.

####################
# Build Dataframes #
####################
def buildDF(inputDir, sampleList):
    if len(sampleList) == 1:
        sample = inputDir + '/' + sampleList[0] + '/' +'total_filtered.vcf'
        dataframe = buildSingleDF(sample)

        return dataframe

    else:
        dataframe = buildAverageStructure(inputDir, sampleList)

        return dataframe

################################
# Single Sample Data Structure #
################################
def buildSingleDF(sample):
    from parseLine import parseLine
    target = open(sample, 'r')
    dataFrame = {}

    for line in target:
        if '#' not in i and 'chr' in i: # skip vcf info
            chrom, loc, AONum, DPNum, var, WT, AFNum = parseLine(line)
            location = '%s-%s-%s' % (chrom, str(loc), str(var))

            dataframe[variant] = {'var': var, 'vaf': AFNum, 'chr': chrom, 'wt': WT}

    return dataframe

############################
# Build Avg Data Structure #
############################
def buildAverageStructure(inputDir, samples):
    from parseLine import parseLine
    tempData = {}
    for i in samples:
        target = open(inputDir + '/' + i + '/' +'total_filtered.vcf', 'r')
        for line in target:
            if '#' not in line and 'chr' in line: # skip the info
                chrom, loc, AONum, DPNum, var, WT, AFNum = parseLine(line)
                # Ex: location = chr1:1234:A
                location = '%s-%s-%s' % (chrom, str(loc), str(var))

                if location in tempData:
                    tempData[location]['vaf'].append(AFNum)
                else:
                    tempData[location] = {'vaf':[AFNum]}

    # average all of the AFNum values
    avgData = takeAverage(tempData)
    return avgData

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
