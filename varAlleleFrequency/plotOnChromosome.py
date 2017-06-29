#!/usr/bin/env python
# the purpose of this script is to transform vafRepeatability.txt
# so it can be plotted in R because R is stoopid

def createInput():
    inFile = open('vafRepeatability.txt', 'r')
    outFile = open('vafPlotting.txt', 'w')

    count = 0
    for line in inFile:
        if count == 0:
            count+=1
            outFile.write('vaf\tchrom\tlocus\twt\tmut\n')
        else:
            vaf=line.split('\t')[1]
            x=line.split('\t')[2].strip('\n').split('-')
            chrom = x[0].strip('chr')
            locus = float(x[1])
            wt = x[2]
            mut = x[3]
            outFile.write('%s\t%s\t%s\t%s\t%s\n' % (vaf, chrom,locus,wt,mut))

    inFile.close()
    outFile.close()

def plot():
    inFile = open('vafPlotting.txt', 'r')
    vafs=[]
    loci=[]
    changes=[]

    # {(chrom,loc):{chrom:4,locus:[1234],vaf:[0.5],mut:[CT]}}
    totalDict = {}

    count = 0
    for line in inFile:
        # skip headers
        if count == 0:
            count += 1
        else:
            x = line.split('\t')
            chrom = x[1]
            loc = float(x[2])
            wt = x[3]
            mut = x[4].strip('\n')
            vaf = float(x[0])

            #if chrom == '4':
            totalDict = checkLoc(chrom, loc, wt, mut, vaf, totalDict)

    for probe in totalDict:
        generatePlot(totalDict[probe])

# this will check if region is in a currently probed region or
# in a differently probed region
def checkLoc(chrom, loc, wt, mut, vaf, totalDict):
    # total Dict form:
    # {(chrom,loc):{chrom:4,locus:[1234],vaf:[0.5],mut:[CT]}}

    key = (chrom, loc)
    isunique=True

    for i in totalDict:
        # is location prob in same probe
        if chrom == i[0] and loc < i[1] + 200 and loc > i[1] - 200:
            key = i
            isunique=False

    if isunique:
        totalDict[key]={'chrom':0,'locus':[],'vaf':[],'mut':[]}
        totalDict[key]['chrom']=chrom
        totalDict[key]['locus']=[loc]
        totalDict[key]['vaf']=[vaf]
        totalDict[key]['mut']=[('%s%s' % (wt, mut))]
    else:
        totalDict[key]['locus'].append(loc)
        totalDict[key]['vaf'].append(vaf)
        totalDict[key]['mut'].append(('%s%s' % (wt, mut)))

    return totalDict

def generateColors(muts):
    from Bio.Seq import Seq
    # same colors as base bias plotting
    colors = {'CA':'c', 'CG':'k', 'CT':'r', 'TA':'0.5', 'TC':'g', 'TG':'m'}
    changes = []

    # change to cononical changes
    for i in muts:
        if i in colors:
            change = i
        else:
            change = str(Seq(i).complement())

        changes.append(colors[change])

    return changes

def identifyLocus(probe):
    ID = ''
    probes = {'TIIIa':['1152278','1152279'],'NRAS-1':['1152564','1152565','1152566'],'NRAS-2':['1152587','1152588'],'DNMT3a':['254572','254573'],'IDH1':['2091130','2091131','2091132'],'SF3B1':['1982668','1982669'],'TIIIb':['2231906','2231907','2231908'],'TIIIc':['2290411','2290412'],'TET2-1':['1061972','1061973','1061974'],'TET2-2':['1061551','1061552'],'TIIId':['1105411','1105412','1105413'],'TIIIe':['1129972','1129973'],'TIIIf':['1211677','1211678'],'TIIIg':['1235477','1235478','1235479'],'TIIIh':['1244286','1244287'],'JAK2':['50737','50738'],'TIIIj':['21262','21263','21264'],'TIIIk':['2389','2390'],'TIIIl':['2593','2594'],'TIIIm':['114865','114866','114867'],'HRAS':['5342','5343'],'KRAS-1':['253982','253983','253984'],'KRAS-2':['253802','253803'],'TIIIn':['925270','925271'],'IDH2':['906318','906319'],'TIIIo':['733796','733797','733798'],'TIIIp':['824550','824551'],'TIIIq':['859491','859492'],'p53-1':['75775','75776'],'p53-2':['75783','75784','75785'],'p53-3':['75770','75771','75772'],'GATA1':['486496','486497','486498']}

    for i in probes:
        for j in probes[i]:
            if j in probe:
                ID = i

    return ID

def generatePlot(probe):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    # probe in form:
    # {locus:[1234],vaf:[0.5],mut:[CT]}
    muts = probe['mut']
    changes = generateColors(muts)
    loci = probe['locus']
    vafs = probe['vaf']

    # create plot
    ind = np.arange(len(loci))
    p1 = plt.bar(loci, vafs, color=changes)
    # when youre bored
    #plt.xkcd()

    # ID gene target
    gene = identifyLocus(str(probe['locus'][0]))

    plt.ylim(ymin=0)
    plt.ylim(ymax=0.003)
    plt.xlabel('Chromosome Location')
    plt.ylabel('Variant Allele Frequencies')
    plt.title('Probe Within Chromosome %s:%s' % (probe['chrom'],gene))

    # this formats the xticks so they aren't abbreviated
    from matplotlib.ticker import FormatStrFormatter
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%d'))
    plt.xticks(rotation=45)

    # create the legend
    cyan = mpatches.Patch(color='cyan', label='C-A')
    black = mpatches.Patch(color='black', label='C-G')
    red = mpatches.Patch(color='red', label='C-T')
    gray = mpatches.Patch(color='gray', label='T-A')
    green = mpatches.Patch(color='green', label='T-C')
    magenta = mpatches.Patch(color='magenta', label='T-G')
    plt.legend(handles=[cyan, black, red, gray, green, magenta])

    plt.show()

if __name__ == '__main__':
    createInput()
    plot()
    #plot('2', 27000000)
    #plot('16', 74379750)
    #plot('4', 134428737)
    #plot('11', 2226328)
