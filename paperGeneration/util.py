#!/usr/bin/env python

def defineProbes():
    probes = {'TIIIa':['1152278','1152279'],'NRAS-1':['1152564','1152565','1152566'],'NRAS-2':['1152587','1152588'],'DNMT3a':['254572','254573'],'IDH1':['2091130','2091131','2091132'],'SF3B1':['1982668','1982669'],'TIIIb':['2231906','2231907','2231908'],'TIIIc':['2290411','2290412'],'TET2-1':['1061972','1061973','1061974'],'TET2-2':['1061551','1061552'],'TIIId':['1105411','1105412','1105413'],'TIIIe':['1129972','1129973'],'TIIIf':['1211677','1211678'],'TIIIg':['1235477','1235478','1235479'],'TIIIh':['1244286','1244287'],'JAK2':['50737','50738'],'TIIIj':['21262','21263','21264'],'TIIIk':['2390'],'TIIIl':['2593','2594'],'TIIIm':['114865','114866','114867'],'HRAS':['5342','5343'],'KRAS-1':['253982','253983','253984'],'KRAS-2':['253802','253803'],'TIIIn':['925270','925271'],'IDH2':['906318','906319'],'TIIIo':['733796','733797','733798'],'TIIIp':['824550','824551'],'TIIIq':['859491','859492'],'p53-1':['75775','75776'],'p53-2':['75783','75784','75785'],'p53-3':['75770','75771','75772'],'GATA1':['486496','486497','486498']}
    return probes

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



# gets sequencing surrounding a particular position using a reference genome
def getRefSequence(vcfObj, upstream, downstream, ref):
    #from parseline import VCFObj
    from subprocess import check_output, STDOUT
    from string import upper

    low = int(vcfObj.location) - upstream
    high = int(vcfObj.location) + downstream
    temp = check_output('samtools faidx %s %s:%s-%s' % (ref, vcfObj.chrom, low, high), stderr=STDOUT, shell=True)

    finalSeq = ''
    for line in temp.split('\n'):
        if '>' not in line:
            finalSeq += line

    finalSeq = finalSeq.upper()
    return finalSeq

# this creates a general pandas DataFrame that can be used for a lot of different analysis
'''
         Loc WT Var Change ConvChange  AO     DP          VAF IntEx  \
0      115227854  T   A    T>A        T>A   5  45143  0.000110759  TIII
1      115227854  T   C    T>C        T>C  45  45143  0.000996832  TIII
2      115227855  T   A    T>A        T>A  15  45199  0.000331866  TIII
3      115227855  T   C    T>C        T>C  32  45199   0.00070798  TIII

Upstream  Downstream Individual
0      GCTTCCCTGC  TCCAATGGCA          1
1      GCTTCCCTGC  TCCAATGGCA          1
2      CTTCCCTGCT  CCAATGGCAT          1
3      CTTCCCTGCT  CCAATGGCAT          1
'''

def populatePandasDataframe(dirinput, fileList, probes, ref, upstream=10, downstream=10):
    import pandas as pd
    from Bio.Seq import Seq
    print('Building data structure...')

    allSamples = []
    columns = ['Chrom','Loc','WT','Var','Change','ConvChange','AO','DP','VAF','IntEx','Gene','Upstream','Downstream','Individual']
    dat = []

    tempAllVariants = []
    sampleCount = 0
    for sample in fileList:
        inFile = open(dirinput + '/' + sample + '/onlyProbedRegions.vcf', 'r')
        sampleCount += 1

        for line in inFile:
            if '#' not in line and 'chr' in line: # skip the info
                lineobj = VCFObj(line)
                # convert to six changes
                if lineobj.wt == 'G' or lineobj.wt == 'A':
                    wt = str(Seq(lineobj.wt).complement())
                    var = str(Seq(lineobj.var).complement())
                else:
                    wt = str(lineobj.wt)
                    var = str(lineobj.var)

                surrounding = getRefSequence(lineobj, upstream, downstream, ref)
                up = str(surrounding[:upstream])
                down = str(surrounding[-downstream:])

                probeRegion = ''
                for probe in probes:
                    if len(probeRegion) < 1:
                        for loc in probes[probe]:
                            if str(loc) in str(lineobj.location):
                                if probe[0] == 'T':
                                    probeRegion = 'TIII'
                                    specificProbe = probe
                                else:
                                    probeRegion = 'Exon'
                                    specificProbe = probe


                if len(lineobj.wt) == 1 and len(lineobj.var) == 1:
                    dat = [lineobj.chrom, lineobj.location, str(lineobj.wt), str(lineobj.var), str(lineobj.wt) + '>' + str(lineobj.var), wt + '>' + var, lineobj.ao, lineobj.dp, lineobj.af, probeRegion, specificProbe, up, down, sampleCount]
                    tempdat = pd.DataFrame(dat, index=columns)
                    tempAllVariants.append(tempdat.T)

        inFile.close()
    allVariants = pd.concat(tempAllVariants, ignore_index=True)

    return allVariants

# multipurpose plotting
def plotData(means, stddev, colors, labels, xlabel, ylabel, title):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    axis_font = {'fontname':'Arial', 'size':'30'}
    tick_font = {'fontname':'Arial', 'size':'30'}
    plt.bar(range(len(means)), means, color=colors)
    plt.errorbar(range(len(stddev)), means, yerr=stddev, linestyle='None', color='black')
    plt.xticks(range(len(labels)), labels, rotation=90, **tick_font)
    plt.yticks(**tick_font)
    plt.ylabel(ylabel, **axis_font)
    plt.xlabel(xlabel, **axis_font)
    plt.title(title, **axis_font)
    #exon = mpatches.Patch(color='blue', label='Exon')
    #intron = mpatches.Patch(color='orange', label='Intron')
    #plt.legend(handles=[exon, intron])
    plt.show()

def plotStacked(means, std, xlabel, ylabel, title):
    import matplotlib.pyplot as plt
    import pandas as pd

    axis_font = {'fontname':'Arial', 'size':'30'}
    tick_font = {'fontname':'Arial', 'size':'30'}
    meandf=pd.DataFrame(means)
    stddf = pd.DataFrame(std)
    positions = [-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10]

    meandf.plot(kind='bar', yerr=stddf, stacked=False)
    plt.xticks(range(len(positions)), positions, **tick_font)
    plt.yticks(**tick_font)
    #plt.ylim(0,0.75)
    plt.ylabel(ylabel, **axis_font)
    plt.xlabel(xlabel, **axis_font)
    plt.title(title, **axis_font)
    plt.show()

    meandf.plot(kind='bar', yerr=stddf, stacked=True)
    plt.xticks(range(len(positions)), positions, **tick_font)
    plt.yticks(**tick_font)
    plt.ylabel(ylabel, **axis_font)
    plt.xlabel(xlabel, **axis_font)
    plt.title(title, **axis_font)
    plt.show()

def saveData(df, name):
    import pickle
    print('Saving data...')
    output = open(name, 'wb')
    pickle.dump(df, output)
    output.close()

def loadData(name):
    import pickle
    print('Loading data...')
    target = open(name, 'rb')
    df = pickle.load(target)
    target.close()
    return df

def overallCaptureFractionPandas(allVariants, samples):
    import numpy as np
    fractions = {'C>A':[], 'C>G':[], 'C>T':[], 'T>A':[], 'T>C':[], 'T>G':[]}
    order = ['C>A','C>G','C>T','T>A','T>C','T>G']
    means = []
    std = []

    for indiv in range(len(samples)):
        totAO = allVariants[allVariants['Individual'] == indiv + 1]['AO'].sum()
        for change in fractions:
            x = allVariants[allVariants['Individual'] == indiv + 1]
            ao = x[x['ConvChange'] == change]['AO'].sum()
            fraction = float(ao) / totAO
            fractions[change].append(fraction)

    for change in order:
        means.append(np.mean(fractions[change]))
        std.append(np.std(fractions[change]))

    return means, std, order

def plotSixChanges(means, stddev, colors, labels, xlabel, ylabel, title, output):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import seaborn as sns

    #colors=['cyan','black','red','grey','green','pink']
    sns.barplot(range(len(means)), means, palette=colors)
    plt.errorbar(range(len(stddev)), means, yerr=stddev, linestyle='None', color='#575859')
    plt.xlabel(xlabel, {'size':'20'})
    plt.ylabel(ylabel, {'size':'20'})
    plt.title(title, {'size':'20'})
    #plt.ylim(0,0.5)

    sns.set_context("paper", font_scale=1.1)
    sns.despine(offset=10, trim=True)
    plt.xticks(range(len(labels)), labels, rotation=90)
    plt.tight_layout()
    plt.savefig(output, bbox_inches='tight',dpi=1200)
    plt.show()

def quantifyCpGVAFs(allVariants, samples):
    import numpy as np
    cpgmeans = []
    nonmeans = []

    mean = []
    std = []

    x = allVariants[allVariants['WT'] == 'C']
    cpg = x[x.Downstream.str[0] == ('G')]
    noncpg = x[x.Downstream.str[0] != ('G')]

    for indiv in range(len(samples)):
        cpgmean = cpg['VAF'].mean()
        nonmean = noncpg['VAF'].mean()

        cpgmeans.append(cpgmean)
        nonmeans.append(nonmean)

    mean.append(np.mean(nonmeans))
    std.append(np.std(nonmeans))
    mean.append(np.mean(cpgmeans))
    std.append(np.std(cpgmeans))

    return mean, std

def overallFractionByTriplet(allVariants, samples):
    import numpy as np
    from collections import defaultdict

    means = []
    std = []
    temp = defaultdict(list)

    changes = ['C>A','C>G','C>T','T>A','T>C','T>G']
    order = ['TCT','CCT','GCT','ACT','TCC','CCC','GCC','ACC','TCG','CCG','GCG','ACG','TCA','CCA','GCA','ACA',
            'TTT','CTT','GTT','ATT','TTC','CTC','GTC','ATC','TTG','CTG','GTG','ATG','TTA','CTA','GTA','ATA']
    labels = [
            'TCT','CCT','GCT','ACT','TCC','CCC','GCC','ACC','TCG','CCG','GCG','ACG','TCA','CCA','GCA','ACA',
            'TCT','CCT','GCT','ACT','TCC','CCC','GCC','ACC','TCG','CCG','GCG','ACG','TCA','CCA','GCA','ACA',
            'TCT','CCT','GCT','ACT','TCC','CCC','GCC','ACC','TCG','CCG','GCG','ACG','TCA','CCA','GCA','ACA',
            'TTT','CTT','GTT','ATT','TTC','CTC','GTC','ATC','TTG','CTG','GTG','ATG','TTA','CTA','GTA','ATA',
            'TTT','CTT','GTT','ATT','TTC','CTC','GTC','ATC','TTG','CTG','GTG','ATG','TTA','CTA','GTA','ATA',
            'TTT','CTT','GTT','ATT','TTC','CTC','GTC','ATC','TTG','CTG','GTG','ATG','TTA','CTA','GTA','ATA',
            ]
    colors=['#4C72B0','#4C72B0','#4C72B0','#4C72B0','#4C72B0','#4C72B0','#4C72B0','#4C72B0','#4C72B0','#4C72B0','#4C72B0','#4C72B0','#4C72B0','#4C72B0','#4C72B0','#4C72B0',
            '#55A868','#55A868','#55A868','#55A868','#55A868','#55A868','#55A868','#55A868','#55A868','#55A868','#55A868','#55A868','#55A868','#55A868','#55A868','#55A868',
            '#C44E52','#C44E52','#C44E52','#C44E52','#C44E52','#C44E52','#C44E52','#C44E52','#C44E52','#C44E52','#C44E52','#C44E52','#C44E52','#C44E52','#C44E52','#C44E52',
            '#8172B2','#8172B2','#8172B2','#8172B2','#8172B2','#8172B2','#8172B2','#8172B2','#8172B2','#8172B2','#8172B2','#8172B2','#8172B2','#8172B2','#8172B2','#8172B2',
            '#CCB974','#CCB974','#CCB974','#CCB974','#CCB974','#CCB974','#CCB974','#CCB974','#CCB974','#CCB974','#CCB974','#CCB974','#CCB974','#CCB974','#CCB974','#CCB974',
            '#64B5CD','#64B5CD','#64B5CD','#64B5CD','#64B5CD','#64B5CD','#64B5CD','#64B5CD','#64B5CD','#64B5CD','#64B5CD','#64B5CD','#64B5CD','#64B5CD','#64B5CD','#64B5CD'
            ]


    for indiv in range(len(samples)):
        aotot = allVariants[allVariants['Individual'] == indiv + 1]['AO'].sum()
        z = allVariants[allVariants['Individual'] == indiv + 1]
        for change in changes:
            a = z[z.ConvChange.str[0] == change[0]]
            b = a[a.ConvChange.str[2] == change[2]]
            for triplet in order:
                c = b[b.Upstream.str[-1] == triplet[0]]
                d = c[c.Downstream.str[0] == triplet[2]]['AO'].sum()

                temp[(change, triplet)].append(float(d) / aotot)

    for change in changes:
        for triplet in order:
            if triplet[1] == change[0]:
                means.append(np.mean(temp[(change, triplet)]))
                std.append(np.std(temp[(change, triplet)]))

    return means, std, labels, changes, colors
























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

