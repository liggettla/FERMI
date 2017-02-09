#!/usr/bin/env python
# The purpose of this script is to compute an R-squared
# value for the vaf repeatability plots

def getCoords(line):
    columns = line.split('\t')
    x = float(columns[0])
    y = float(columns[1])
    loc = columns[2].strip('\n')

    return x, y, loc

# The assumption is that x == y such that predicted values
# should be identical
def computeError(x, y):
    error = y - x
    sqrError = pow(error, 2)

    return sqrError

def computeMean(parsedFile):
    total = 0
    count = 0

    for i in parsedFile:
        count += 1
        total += parsedFile[i]['y']

    yMean = total / count
    return yMean

def meanDistance(yMean, y):
    meanDist = y - yMean
    sqrMeanDist = pow(meanDist, 2)

    return sqrMeanDist

def computeRSquared(parsedFile):
    errorSum = 0
    meanSum = 0

    for i in parsedFile:
        errorSum = parsedFile[i]['sqrError']
        meanSum = parsedFile[i]['sqrMeanDist']

    r2 = 1 - errorSum / meanSum
    return r2

def getRSquared(target='vafRepeatability.txt'):
    from math import pow

    target = 'vafRepeatability.txt'
    inputFile = open(target, 'r')
    parsedFile = {}

    count = 0
    for line in inputFile: # skip headers
        if count == 0:
            count += 1
        else:
            x, y, loc = getCoords(line)
            sqrError = computeError(x, y)
            parsedFile[loc]={'x':x, 'y':y, 'sqrError':sqrError}

    yMean = computeMean(parsedFile)
    for i in parsedFile:
        y = parsedFile[i]['y']
        sqrMeanDist = meanDistance(yMean, y)
        parsedFile[i]['sqrMeanDist'] = sqrMeanDist
    r2 = computeRSquared(parsedFile)
    print r2

if __name__ == '__main__':
    getRSquared()
