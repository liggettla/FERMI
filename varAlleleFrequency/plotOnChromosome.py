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
    import numpy as np
    import matplotlib.pyplot as plt
    from Bio.Seq import Seq
    import matplotlib.patches as mpatches

    inFile = open('vafPlotting.txt', 'r')
    vafs=[]
    loci=[]
    changes=[]

    # form: {(chrom,loc):{vaf:[],mut[],wt[]}}
    totalDict = {}

    # same colors as base bias plotting
    colors = {'CA':'c', 'CG':'k', 'CT':'r', 'TA':'0.5', 'TC':'g', 'TG':'m'}

    for line in inFile:

        x = line.split('\t')
        if x[1] == chrom and float(x[2]) < loc:
            vafs.append(float(x[0]))
            loci.append(float(x[2]))
            mut = x[4].strip('\n')
            wt = x[3]

            if ('%s%s' % (wt, mut)) in colors:
                change = ('%s%s' % (wt, mut))
            else:
                change = str(Seq(('%s%s' % (wt, mut))).complement())

            changes.append(colors[change])

    # create plot
    ind = np.arange(len(loci))
    p1 = plt.bar(loci, vafs, color=changes)
    plt.ylim(ymin=0)
    plt.ylim(ymax=0.003)
    plt.xlabel('Chromosome Location')
    plt.ylabel('Variant Allele Frequencies')

    # create the legend
    cyan = mpatches.Patch(color='cyan', label='C-A')
    black = mpatches.Patch(color='black', label='C-G')
    red = mpatches.Patch(color='red', label='C-T')
    gray = mpatches.Patch(color='gray', label='T-A')
    green = mpatches.Patch(color='green', label='T-G')
    magenta = mpatches.Patch(color='magenta', label='T-G')
    plt.legend(handles=[cyan, black, red, gray, green, magenta])

    plt.show()

# this will check if region is in a currently probed region or
# in a differently probed region
def checkLoc(chrom, loc, wt, mut, vaf, totalDict):
    # total Dict form:
    # {(chrom,loc):{locus:[1234],vaf:[0.5],mut[CT]}}

    key = (chrom, loc)
    isunique=True

    for i in totalDict:
        # is location prob in same probe
        if chrom == i[0] and loc < i[1] + 200 and loc > i[1] - 200:
            key = i
            isunique=False

    if isunique:
        totalDict[key]['locus']=[locus]
        totalDict[key]['vaf']=[vaf]
        totalDict[key]['mut']=[('%s%s' % (wt, mut))]
    else:
        totalDict[key]['locus'].append(locus)
        totalDict[key]['vaf'].append(vaf)
        totalDict[key]['mut'].append(('%s%s' % (wt, mut)))

    return totalDict

if __name__ == '__main__':
    createInput()
    plot('2', 27000000)
    #plot('16', 74379750)
    #plot('4', 134428737)
    #plot('11', 2226328)
