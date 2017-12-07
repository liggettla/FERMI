#!/usr/bin/env python

# this computs r2 and pearson correlation coeff
def computersquared(parsedFile):
    from scipy import stats
    import numpy as np
    from scipy.stats.stats import pearsonr
    x = parsedFile['x']
    y = parsedFile['y']

    # p value here is two-sided p-value for a hypothesis test
    # whose null hypothesis is that the slope is zero.
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

    # calculate pearson correlation coefficient
    p = pearsonr(x, y)

    return r_value**2, p

def getCoords(line):
    columns = line.split('\t')
    x = float(columns[0])
    y = float(columns[1])
    loc = columns[2].strip('\n')

    return x, y, loc

def getRSquared(target='vafRepeatability.txt'):
    inputFile = open(target, 'r')
    parsedFile = {'x':[],'y':[],'loc':[]}

    count = 0
    for line in inputFile: # skip headers
        if count == 0:
            count += 1
        else:
            x, y, loc = getCoords(line)
            parsedFile['x'].append(x)
            parsedFile['y'].append(y)
            parsedFile['loc'].append(loc)

    r2 = computersquared(parsedFile)
    return r2

if __name__ == '__main__':
    r2 = getRSquared()
    print r2

