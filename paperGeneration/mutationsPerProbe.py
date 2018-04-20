#!/usr/bin/env python
# The purpose of this script is to ask how many variants are found within
# a given probed region to see if there is any bias among the probes.

def runArgParse():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--inFile', '-i', type=str, required=True, help='The input vcf to be analyzed.')
    parser.add_argument('--outputDir', '-o', type=str, required=True, help='Output directory for plot file.')

    args = parser.parse_args()
    inFile = args.inFile
    outputDir = args.outputDir

    return inFile, outputDir

def identifyProbe(location):
    probes = {'TIIIa':['1152278','1152279'],'NRAS-1':['1152564','1152565','1152566'],'NRAS-2':['1152587','1152588'],'DNMT3a':['254572','254573'],'IDH1':['2091130','2091131','2091132'],'SF3B1':['1982668','1982669'],'TIIIb':['2231906','2231907','2231908'],'TIIIc':['2290411','2290412'],'TET2-1':['1061972','1061973','1061974'],'TET2-2':['1061551','1061552'],'TIIId':['1105411','1105412','1105413'],'TIIIe':['1129972','1129973'],'TIIIf':['1211677','1211678'],'TIIIg':['1235477','1235478','1235479'],'TIIIh':['1244286','1244287'],'JAK2':['50737','50738'],'TIIIj':['21262','21263','21264'],'TIIIk':['2390'],'TIIIl':['2593','2594'],'TIIIm':['114865','114866','114867'],'HRAS':['5342','5343'],'KRAS-1':['253982','253983','253984'],'KRAS-2':['253802','253803'],'TIIIn':['925270','925271'],'IDH2':['906318','906319'],'TIIIo':['733796','733797','733798'],'TIIIp':['824550','824551'],'TIIIq':['859491','859492'],'p53-1':['75775','75776'],'p53-2':['75783','75784','75785'],'p53-3':['75770','75771','75772'],'GATA1':['486496','486497','486498']}

    for i in probes:
        for j in probes[i]:
            if j in location:
                return i
    return False


def displayTally(tally):
    print('Probe\tNum Vars')
    for i in tally:
        print('%s\t%s' % (i, tally[i]))

def plotTally(outputDir, tally, pdf, title='Variants Per Probe', xlabel='Number of Variants'):
    import matplotlib.pyplot as plt

    plt.bar(range(len(tally)), tally.values(), align='center')
    plt.title(title)
    plt.xlabel('Probes', fontsize=1)
    plt.ylabel(xlabel)
    plt.xticks(range(len(tally)), tally.keys(), rotation='90')
    plt.tight_layout()


    #fig = plt.figure()
    #plt.savefig(output, format='png')
    pdf.savefig()
    #plt.show()

    return pdf

# Averages each of the depths for all variants at a probed region
def getMeanCoverage(totalCoverage):
    from numpy import mean
    for i in totalCoverage:
        totalCoverage[i] = mean(totalCoverage[i])

    return totalCoverage

# This normalizes the number of variants in a probed region my the
# average number of unique captures of that region
def normalizeCounts(raw, totalCoverage):
    for i in raw:
        raw[i] = raw[i] / totalCoverage[i]

    return raw

def mutationsPerProbe(inFile, outputDir):
    from parseLine import seqRead
    target = open(inFile, 'r')

    # Number of unique variants found within a particular capture region
    uniqVars = {'TIIIa':0,'NRAS-1':0,'NRAS-2':0,'DNMT3a':0,'IDH1':0,'SF3B1':0,'TIIIb':0,'TIIIc':0,'TET2-1':0,'TET2-2':0,'TIIId':0,'TIIIe':0,'TIIIf':0,'TIIIg':0,'TIIIh':0,'JAK2':0,'TIIIj':0,'TIIIk':0,'TIIIl':0,'TIIIm':0,'HRAS':0,'KRAS-1':0,'KRAS-2':0,'TIIIn':0,'IDH2':0,'TIIIo':0,'TIIIp':0,'TIIIq':0,'p53-1':0,'p53-2':0,'p53-3':0,'GATA1':0}
    # Number of total variants found within a particular capture region
    totalVars = {'TIIIa':0,'NRAS-1':0,'NRAS-2':0,'DNMT3a':0,'IDH1':0,'SF3B1':0,'TIIIb':0,'TIIIc':0,'TET2-1':0,'TET2-2':0,'TIIId':0,'TIIIe':0,'TIIIf':0,'TIIIg':0,'TIIIh':0,'JAK2':0,'TIIIj':0,'TIIIk':0,'TIIIl':0,'TIIIm':0,'HRAS':0,'KRAS-1':0,'KRAS-2':0,'TIIIn':0,'IDH2':0,'TIIIo':0,'TIIIp':0,'TIIIq':0,'p53-1':0,'p53-2':0,'p53-3':0,'GATA1':0}
    # Number of total probes capturing a particular region
    totalCoverage = {'TIIIa':[],'NRAS-1':[],'NRAS-2':[],'DNMT3a':[],'IDH1':[],'SF3B1':[],'TIIIb':[],'TIIIc':[],'TET2-1':[],'TET2-2':[],'TIIId':[],'TIIIe':[],'TIIIf':[],'TIIIg':[],'TIIIh':[],'JAK2':[],'TIIIj':[],'TIIIk':[],'TIIIl':[],'TIIIm':[],'HRAS':[],'KRAS-1':[],'KRAS-2':[],'TIIIn':[],'IDH2':[],'TIIIo':[],'TIIIp':[],'TIIIq':[],'p53-1':[],'p53-2':[],'p53-3':[],'GATA1':[]}

    for line in target:
        if '#' not in line and 'chr' in line:
            lineObj = seqRead(line)
            probeNum = identifyProbe(lineObj.loc())
            if probeNum:
                uniqVars[probeNum] += 1
                if lineObj.af() < 0.4: # Don't want germline bias
                    totalVars[probeNum] += lineObj.ao()
                totalCoverage[probeNum].append(lineObj.dp())

    print totalCoverage
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(outputDir + '/probeBias.pdf')

    totalCoverage = getMeanCoverage(totalCoverage)
    pdf = plotTally(outputDir, totalCoverage, pdf, 'Probe Bias', 'Avg Num of Captures')

    uniqVars = normalizeCounts(uniqVars, totalCoverage)
    #displayTally(uniqVars)
    pdf = plotTally(outputDir, uniqVars, pdf, 'Normalized Unique Variants Per Probe', 'Number of Variants')

    totalVars = normalizeCounts(totalVars, totalCoverage)
    #displayTally(totalVars)
    pdf = plotTally(outputDir, totalVars, pdf, 'Normalized Total Variants Per Probe', 'Number of Variants')

    pdf.close()

if __name__ == '__main__':
    from parseLine import seqRead
    inFile, outputDir = runArgParse()
    mutationsPerProbe(inFile, outputDir)
