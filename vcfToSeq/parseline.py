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

# this takes a vcf line and outputs a parsed object
class VCFObj:

    def __init__(self, vcfLine):
        self.chrom = str(vcfLine.split('\t')[0])
        self.location = int(str(vcfLine.split('\t')[1]))

        AO = vcfLine.split(';')[5]
        DP = vcfLine.split(';')[7]
        AONum = float(AO.split(',')[0][3:])
        DPNum = float(DP.split(',')[0][3:])
        AFNum = AONum / DPNum

        self.ao = AONum
        self.dp = DPNum
        self.af = AFNum

        self.wt = vcfLine.split('\t')[3]
        self.var = vcfLine.split('\t')[4] # the observed bp change


