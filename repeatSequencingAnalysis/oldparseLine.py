# the purpose of this script is to provide vcf parsing functionality

##############
# Parse Line #
##############
def parseLine(line):
    chrom = str(line.split('\t')[0])
    loc = str(line.split('\t')[1])
    AO = line.split(';')[5]
    DP = line.split(';')[7]
    var = line.split('\t')[4] # the observed bp change
    WT = line.split('\t')[3]

    AONum = float(AO.split(',')[0][3:])
    DPNum = float(DP.split(',')[0][3:])
    AFNum = AONum / DPNum

    return chrom, loc, AONum, DPNum, var, WT, AFNum
