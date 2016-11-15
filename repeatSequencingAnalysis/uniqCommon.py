#!/usr/bin/python
# The purpose of this script is to read in three different frames: a population
# average, sample1 and sample2, where samples 1/2 are technical replicates, in that
# they are independent preparations of the same bleed.
# This script then returns the differences and similarities between these three dataframes.
from itertools import combinations

def uniqCommon(dfAvg, df1, df2):
    dfAvgComp = {'AvgU':{}, 'AvgC':{}, 'RSeqU':{}, 'RSeqC':{}, 'AvgURSeqC':{}, 'AvgURSeqU':{}}
    df1Comp = {'AvgU':{}, 'AvgC':{}, 'RSeqU':{}, 'RSeqC':{}, 'AvgURSeqC':{}, 'AvgURSeqU':{}}
    df2Comp = {'AvgU':{}, 'AvgC':{}, 'RSeqU':{}, 'RSeqC':{}, 'AvgURSeqC':{}, 'AvgURSeqU':{}}

    # below is the structure of input dataframes
    # '{chr16-82455094-C-A': {'var': 'A', 'loc': '82455094', 'chr': 'chr16', 'wt': 'C', 'vaf': 0.00073987950533770217}}

    # structure of output dataframes
    # {AvgU:AvgU[], AvgC:AvgC[], RSeqU:ReSeqU[], RSeqC:ReSeqC[], AvgURSeqC:AvgURSeqC[], AvgURSeqU:AvgURSeqU[]}

    dfList = [df1, df2]
    for a, b in combinations(dfList, 2):
        # associate the dComps with the right incoming df
        if a == df1:
            aComp = df1Comp
            bComp = df2Comp
        elif a == df2:
            aComp = df2Comp
            bComp = df1Comp

        buildDF(a, aComp, b, bComp, dfAvg)

    #dfAvgComp is currently unused
    return dfAvgComp, df1Comp, df2Comp

# this populates both df1Comp and df2Comp
def buildDF(a, aComp, b, bComp, dfAvg):
    # in first iteration a is df1, aComp is df1Comp
    # and b is df2, bComp is df2Comp
    for i in a:
        # compare to dfAvg
        if compareDF(i, dfAvg):
            aComp['AvgC'][i]=a[i]
        else:
            aComp['AvgU'][i]=a[i]

        # compare samples to each other
        if compareDF(i, b):
            aComp['RSeqC'][i]=a[i]
        else:
            aComp['RSeqU'][i]=a[i]

        # compare dfAvg uniques to other sample
        for i in aComp['AvgU']:
            if compareDF(i, b):
                aComp['AvgURSeqC'][i]=a[i]
            else:
                aComp['AvgURSeqU'][i]=a[i]

# this just asks if variant is in common or unique
def compareDF(i, df):
    if i in df:
        return True
    else:
        return False
