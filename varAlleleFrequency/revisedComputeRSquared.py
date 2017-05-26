#!/usr/bin/env python

def computersquared(parsedFile):
    from scipy import stats
    import numpy as np
    x = parsedFile['x']
    y = parsedFile['y']

    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

    return r_value**2

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

