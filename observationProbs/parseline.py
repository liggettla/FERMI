#!/usr/bin/env python
# this provides line parsing functionality that is commonly used
def parseLine(i):
    if '#' not in i and 'chr' in i: # skip the info
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
        return location, AFNum, WT, var, loc, chrom

    else:
        return False, False, False, False, False, False

