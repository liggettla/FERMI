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
    from parseLine import seqRead
    target = open(sample, 'r')
    dataframe = {}

    for line in target:
        if '#' not in line and 'chr' in line: # skip vcf info
            lineObj = seqRead(line)
            # Ex: variant = chr1-1234-G-A
            variant = '%s-%s-%s-%s' % (lineObj.chrom(), lineObj.loc(), lineObj.wt(), lineObj.var())

            # ignore snps in the analysis
            if lineObj.af() < 0.4:
                dataframe[variant] = {'var': lineObj.var(), 'vaf': lineObj.af(), 'chr': lineObj.chrom(), 'wt': lineObj.wt(), 'loc': lineObj.loc()}
            else:
                pass

    return dataframe

############################
# Build Avg Data Structure #
############################
def buildAverageStructure(inputDir, samples):
    from parseLine import seqRead
    tempData = {}
    for i in samples:
        target = open(inputDir + '/' + i + '/' +'total_filtered.vcf', 'r')
        for line in target:
            if '#' not in line and 'chr' in line: # skip the info
                lineObj = seqRead(line)
                # Ex: variant = chr1-1234-G-A
                variant = '%s-%s-%s-%s' % (lineObj.chrom(), lineObj.loc(), lineObj.wt(), lineObj.var())

                if variant in tempData:
                    tempData[variant]['vaf'].append(lineObj.af())
                else:
                    tempData[variant] = {'vaf':[lineObj.af()], 'var': lineObj.var(), 'chr': lineObj.chrom(), 'wt': lineObj.wt(), 'loc': lineObj.loc()}

    # average all of the AFNum values
    avgData = takeAverage(tempData)
    return avgData

###################
# Find Avg AFNums #
###################
# computes the average AFNum for each unique variant
# and returns averaged data structure
def takeAverage(tempData):
    from numpy import mean
    avgData = {}
    for loc in tempData:
        x = mean(tempData[loc]['vaf'])
        tempData[loc]['vaf'] = x

    return tempData
