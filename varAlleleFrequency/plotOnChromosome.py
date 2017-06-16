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

def plot(chrom, loc):
    import numpy as np
    import matplotlib.pyplot as plt
    inFile = open('vafPlotting.txt', 'r')
    vafs=[]
    loci=[]
    muts=[]
    wts=[]
    colors = {'C':'b', 'T':'r', 'G':'g', 'A':'y'}

    for line in inFile:

        x = line.split('\t')
        if x[1] == chrom and float(x[2]) < loc:
            vafs.append(float(x[0]))
            loci.append(float(x[2]))
            muts.append(x[4].strip('\n'))
            wts.append(colors[x[3]])

    print loci, muts, wts
    ind = np.arange(len(loci))
    p1 =plt.bar(loci, vafs, color=wts)
    #plt.colorbar(muts)
    plt.ylim(ymin=0)
    plt.ylim(ymax=0.003)
    #plt.axis([0])
    plt.show()





if __name__ == '__main__':
    createInput()
    #plot('2', 27000000)
    plot('16', 74379750)
    #plot('4', 134428737)
    #plot('11', 2226328)
