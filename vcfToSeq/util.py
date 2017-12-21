#!/usr/bin/env python

def defineProbes():
    probes = {'TIIIa':['1152278','1152279'],'NRAS-1':['1152564','1152565','1152566'],'NRAS-2':['1152587','1152588'],'DNMT3a':['254572','254573'],'IDH1':['2091130','2091131','2091132'],'SF3B1':['1982668','1982669'],'TIIIb':['2231906','2231907','2231908'],'TIIIc':['2290411','2290412'],'TET2-1':['1061972','1061973','1061974'],'TET2-2':['1061551','1061552'],'TIIId':['1105411','1105412','1105413'],'TIIIe':['1129972','1129973'],'TIIIf':['1211677','1211678'],'TIIIg':['1235477','1235478','1235479'],'TIIIh':['1244286','1244287'],'JAK2':['50737','50738'],'TIIIj':['21262','21263','21264'],'TIIIk':['2390'],'TIIIl':['2593','2594'],'TIIIm':['114865','114866','114867'],'HRAS':['5342','5343'],'KRAS-1':['253982','253983','253984'],'KRAS-2':['253802','253803'],'TIIIn':['925270','925271'],'IDH2':['906318','906319'],'TIIIo':['733796','733797','733798'],'TIIIp':['824550','824551'],'TIIIq':['859491','859492'],'p53-1':['75775','75776'],'p53-2':['75783','75784','75785'],'p53-3':['75770','75771','75772'],'GATA1':['486496','486497','486498']}
    return probes

# this returns a data structure with fraction of changes by type of change
# this is independent of trinucleotide context
def overallCaptureFraction(dirinput, fileList, probes, region='intron'):
    print('Reading overall captured base changes...')
    from parseline import VCFObj
    from collections import defaultdict
    from Bio.Seq import Seq

    vafs = {'C>A':[], 'C>G':[], 'C>T':[], 'T>A':[], 'T>C':[], 'T>G':[]}

    for i in fileList:
        temp = {'C>A':0, 'C>G':0, 'C>T':0, 'T>A':0, 'T>C':0, 'T>G':0}
        target = open(dirinput + '/' + i, 'r')

        for line in target:
            if '#' not in line and 'chr' in line: # skip the info
                vcfObj = VCFObj(line)

                # only use single substitutions and eliminate SNPs
                if len(vcfObj.wt) == 1 and len(vcfObj.var) == 1 and vcfObj.af < 0.1:
                    varType = ('%s>%s' % (vcfObj.wt, vcfObj.var))

                    if varType not in vafs:
                        varType = ('%s>%s' % (str(Seq(vcfObj.wt).complement()), str(Seq(vcfObj.var).complement())))

                    for i in probes:
                        for j in probes[i]:
                            if j in str(vcfObj.location):
                                if region == 'intron' and i[0] == 'T':
                                    # add AO together
                                    temp[varType] += int(vcfObj.ao)
                                elif region == 'exon' and not i[0] == 'T':
                                    # add AO together
                                    temp[varType] += int(vcfObj.ao)

        target.close()

        # calculate fractions
        totalAO = 0
        for i in temp:
            totalAO += temp[i]

        for i in temp:
            temp[i] = float(temp[i]) / float(totalAO)

        # transfer data into vafs
        for i in temp:
            vafs[i].append(temp[i])

    return vafs

# generates the 6 change plots
def sixChangeFractions(fractions, title='Exons'):
    print('Plotting substitution fractions...')
    import matplotlib.pyplot as plt
    import numpy as np

    order = ['C>A','C>G','C>T','T>A','T>C','T>G']
    ordered = []
    for i in order:
        ordered.append(fractions[i])

    nums = np.arange(len(fractions))
    vafmean = []
    vafstddev = []

    for i in order:
        vafmean.append(np.mean(fractions[i]))
        vafstddev.append(np.std(fractions[i]))

    plt.clf()
    axis_font = {'fontname':'Arial', 'size':'30'}
    colorlist = ['cyan','black','red','grey','green','magenta']
    plt.errorbar(nums, vafmean, yerr=vafstddev, linestyle="None", color='black')
    plt.bar(nums, vafmean, color=colorlist)
    plt.xticks(range(len(fractions)), order, rotation=0, **axis_font)
    plt.yticks(**axis_font)
    plt.ylabel('Fraction of Observations', **axis_font)
    plt.xlabel('Substitution Type', **axis_font)
    plt.title(title, **axis_font)

    plt.tick_params(
    axis='both',       # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    left='off',        # y-axis ticks are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off

    plt.show()

# this will generate the six change plot with exon and intron
# regions plotted side by side
def plotSixExonIntronTogether(intron, exon):
    import matplotlib.pyplot as plt
    import numpy as np

    order = ['C>A','C>G','C>T','T>A','T>C','T>G']
    labels = ['C>A TIII','C>A Exon','C>G TIII','C>G Exon','C>T TIII','C>T Exon','T>A TIII','T>A Exon','T>C TIII','T>C Exon','T>G TIII','T>G Exon']
    means = []
    std = []
    for i in order:
        means.append(np.mean(intron[i]))
        std.append(np.std(intron[i]))
        means.append(np.mean(exon[i]))
        std.append(np.std(exon[i]))

    plt.clf()
    axis_font = {'fontname':'Arial', 'size':'30'}
    colorlist = ['cyan','cyan','black','black','red','red','grey','grey','green','green','magenta','magenta']
    plt.errorbar(range(len(means)), means, yerr=std, linestyle="None", color='black')
    plt.bar(range(len(means)), means, color=colorlist)
    plt.xticks(range(len(means)), labels, rotation=90, **axis_font)
    plt.yticks(**axis_font)
    plt.ylabel('Fraction of Observations', **axis_font)
    plt.xlabel('Substitution Type', **axis_font)
    plt.title('Intron vs Exon', **axis_font)

    plt.show()



