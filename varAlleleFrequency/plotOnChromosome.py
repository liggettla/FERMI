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
    # if combining
    colors = {'CA':'c', 'CG':'k', 'CT':'#ed0000', 'TA':'0.5', 'TC':'g', 'TG':'m'}
    # if not combining
    colors = {'CA':'c', 'CG':'k', 'CT':'#ed0000', 'TA':'0.5', 'TC':'g', 'TG':'m', 'GT':'#000099', 'GC':'#660066', 'GA':'#00ff00', 'AT':'#ffff00', 'AG':'#ff6600', 'AC':'#663300'}
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

    title_font = {'fontname':'Arial', 'size':'16', 'color':'black', 'weight':'normal'}
    axis_font = {'fontname':'Arial', 'size':'30'}

    plt.ylim(ymin=0)
    plt.ylim(ymax=0.003)
    plt.xlabel('Chromosome Location', **axis_font)
    plt.ylabel('Variant Allele Frequencies', **axis_font)
    plt.title('Probe Within Chromosome %s:%s' % (probe['chrom'],gene), **title_font)

    # this formats the xticks so they aren't abbreviated
    from matplotlib.ticker import FormatStrFormatter
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%d'))
    plt.xticks(rotation=45)

    # create the legend
    cyan = mpatches.Patch(color='cyan', label='C-A')
    black = mpatches.Patch(color='black', label='C-G')
    red = mpatches.Patch(color='#ed0000', label='C-T')
    gray = mpatches.Patch(color='gray', label='T-A')
    green = mpatches.Patch(color='green', label='T-C')
    magenta = mpatches.Patch(color='magenta', label='T-G')
    a = mpatches.Patch(color='#000099', label='G-T')
    b = mpatches.Patch(color='#660066', label='G-C')
    c = mpatches.Patch(color='#00ff00', label='G-A')
    d = mpatches.Patch(color='#ffff00', label='A-T')
    e = mpatches.Patch(color='#ff6600', label='A-G')
    f = mpatches.Patch(color='#663300', label='A-C')
    #plt.legend(handles=[cyan, black, red, gray, green, magenta])
    plt.legend(handles=[cyan, black, red, gray, green, magenta, a,b,c,d,e,f])

    plt.show()

def plot():
    inFile = open('vafPlotting.txt', 'r')
    vafs=[]
    loci=[]
    changes=[]

    # {(chrom,):{chrom:4,locus:[1234],vaf:[0.5],mut:[CT]}}
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

def idExactLocus(chrom, loc, wt, mut, vaf, totalDict):
    exactProbes = {'TIIIa':['115227814','115227978'],'NRAS-1':['115256496','115256680'],'NRAS-2':['115258713','115258897'],'DNMT3a':['25457211','25457364'],'IDH1':['209113077','209113239'],'SF3B1':['198266803','198266967'],'TIIIb':['223190674','223190820'],'TIIIc':['229041101','229041289'],'TET2-1':['106197237','106197405'],'TET2-2':['106155137','106155275'],'TIIId':['110541172','110541302'],'TIIIe':['112997214','112997386'],'TIIIf':['121167756','121167884'],'TIIIg':['123547743','123547901'],'TIIIh':['124428637','124428767'],'JAK2':['5073733','5073887'],'TIIIj':['2126256','2126420'],'TIIIk':['2389983','2390171'],'TIIIl':['2593889','2594074'],'TIIIm':['11486596','11486728'],'HRAS':['534258','534385'],'KRAS-1':['25398247','25398415'],'KRAS-2':['25380242','25380368'],'TIIIn':['92527052','92527176'],'IDH2':['90631809','90631969'],'TIIIo':['73379656','73379832'],'TIIIp':['82455026','82455164'],'TIIIq':['85949137','85949299'],'p53-1':['7577504','7577635'],'p53-2':['7578369','7578544'],'p53-3':['7577084','7577214'],'GATA1':['48649667','48649849']}

    # total Dict form:
    # {(probe):{chrom:4,locus:[1234],vaf:[0.5],mut:[CT]}}

    key = ''
    isunique = True

    for probe in exactProbes:
        match = int(loc)-20 > int(exactProbes[probe][0]) and int(loc) < int(exactProbes[probe][1])
        if match:
            key = probe

    if key in totalDict:
        isunique = False

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

def plotExact():
    inFile = open('vafPlotting.txt', 'r')
    vafs=[]
    loci=[]
    changes=[]

    # {(probe):{chrom:4,locus:[1234],vaf:[0.5],mut:[CT]}}
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
            totalDict = idExactLocus(chrom, loc, wt, mut, vaf, totalDict)

    for probe in totalDict:
        exactPlotting(probe, totalDict[probe])

def exactPlotting(gene, probe):
    exactProbes = {'TIIIa':['115227814','115227978'],'NRAS-1':['115256496','115256680'],'NRAS-2':['115258713','115258897'],'DNMT3a':['25457211','25457364'],'IDH1':['209113077','209113239'],'SF3B1':['198266803','198266967'],'TIIIb':['223190674','223190820'],'TIIIc':['229041101','229041289'],'TET2-1':['106197237','106197405'],'TET2-2':['106155137','106155275'],'TIIId':['110541172','110541302'],'TIIIe':['112997214','112997386'],'TIIIf':['121167756','121167884'],'TIIIg':['123547743','123547901'],'TIIIh':['124428637','124428767'],'JAK2':['5073733','5073887'],'TIIIj':['2126256','2126420'],'TIIIk':['2389983','2390171'],'TIIIl':['2593889','2594074'],'TIIIm':['11486596','11486728'],'HRAS':['534258','534385'],'KRAS-1':['25398247','25398415'],'KRAS-2':['25380242','25380368'],'TIIIn':['92527052','92527176'],'IDH2':['90631809','90631969'],'TIIIo':['73379656','73379832'],'TIIIp':['82455026','82455164'],'TIIIq':['85949137','85949299'],'p53-1':['7577504','7577635'],'p53-2':['7578369','7578544'],'p53-3':['7577084','7577214'],'GATA1':['48649667','48649849']}
    # gene = NRAS-1
    # probe = {chrom:4,locus:[1234],vaf:[0.5],mut:[CT]}
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    # probe in form:
    # {locus:[1234],vaf:[0.5],mut:[CT]}
    #print probe['vaf']
    muts = probe['mut']
    changes = generateColors(muts)
    loci = probe['locus']
    vafs = probe['vaf']

    # create plot
    ind = np.arange(len(loci))
    p1 = plt.bar(loci, vafs, color=changes)
    # when youre bored
    #plt.xkcd()

    plt.ylim(ymin=0)
    plt.ylim(ymax=0.003)

    title_font = {'fontname':'Arial', 'size':'30', 'color':'black', 'weight':'normal'}
    axis_font = {'fontname':'Arial', 'size':'30'}

    # there is some weird error were the second element in exactProbes
    # is not a gene name, this bypasses that problem, but better to fix it
    if len(gene) > 3:
        plt.xlim(xmin=int(exactProbes[gene][0]))
        plt.xlim(xmax=int(exactProbes[gene][1]))

    plt.xlabel('Chromosome Location', **axis_font)
    plt.ylabel('Variant Allele Frequencies', **axis_font)
    plt.title('Probe Within Chromosome %s:%s' % (probe['chrom'],gene), **title_font)
    plt.savefig('plotOnProbe.png')

    # this formats the xticks so they aren't abbreviated
    from matplotlib.ticker import FormatStrFormatter
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%d'))
    plt.xticks(rotation=45, **axis_font)
    plt.yticks(**axis_font)


    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off


    # create the legend
    cyan = mpatches.Patch(color='cyan', label='C-A')
    black = mpatches.Patch(color='black', label='C-G')
    red = mpatches.Patch(color='#ed0000', label='C-T')
    gray = mpatches.Patch(color='gray', label='T-A')
    green = mpatches.Patch(color='green', label='T-C')
    magenta = mpatches.Patch(color='magenta', label='T-G')
    a = mpatches.Patch(color='#000099', label='G-T')
    b = mpatches.Patch(color='#660066', label='G-C')
    c = mpatches.Patch(color='#00ff00', label='G-A')
    d = mpatches.Patch(color='#ffff00', label='A-T')
    e = mpatches.Patch(color='#ff6600', label='A-G')
    f = mpatches.Patch(color='#663300', label='A-C')
    #plt.legend(handles=[cyan, black, red, gray, green, magenta])
    plt.legend(handles=[cyan, black, red, gray, green, magenta, a,b,c,d,e,f], prop={'size': 15})
    plt.tight_layout()

    plt.show()

if __name__ == '__main__':
    createInput()
    #plot()
    plotExact()
