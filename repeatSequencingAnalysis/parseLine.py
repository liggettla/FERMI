#!/usr/bin/python
# The purpose of this script is to provide vcf parsing functionality

# This class takes an individual line from a vcf file and outputs important info
class seqRead(object):

    def __init__(self, line):
        self.line = line

    def chrom(self):
        chrom = str(self.line.split('\t')[0])
        return chrom

    def loc(self):
        loc = str(self.line.split('\t')[1])
        return loc

    def wt(self):
        WT = self.line.split('\t')[3]
        return WT

    def var(self):
        var = self.line.split('\t')[4] # the observed bp change
        return var

    def ao(self):
        AO = self.line.split(';')[5]
        AONum = float(AO.split(',')[0][3:])
        return AONum

    def dp(self):
        DP = self.line.split(';')[7]
        DPNum = float(DP.split(',')[0][3:])
        return DPNum

    def af(self):
        AO = self.line.split(';')[5]
        AONum = float(AO.split(',')[0][3:])
        DP = self.line.split(';')[7]
        DPNum = float(DP.split(',')[0][3:])
        AFNum = AONum / DPNum
        return AFNum
