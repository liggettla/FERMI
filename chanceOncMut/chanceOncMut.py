#!/usr/bin/env python

#This purpose of this script is to ask the question of whether or not
#oncogenic mutations are being seen more often chance would predict based on
#the frequency of mutations in TIII and non-oncogenic regions as compared
#with the frequency seen at oncogenic locations
import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--indir', '-i', required=True, type=str, help='Specifies the input directory containing all folders containing output analysis from a fermi analysis run.')
parser.add_argument('--outdir', '-o', required=True, type=str, help='Specifies the output directory location for the analysis output file')
parser.add_argument('--samples', '-s', required=True, nargs='*', type=str, help='Name of the directories containing fermi analysis of samples to be compared.')

args = parser.parse_args()

inputDir = args.indir
outputDir = args.outdir
samples = args.samples # this is a list

oncoSites = [5073770, 7577539,7577119,115256529,115258747,115258744,534287,534288,534289,25398284,25380275,106197266,106197267,106197268,106197269,106155172,106155173,106155174,25457242,25457243,209113112,209113113,90631934,90631838,48649700]
oncoGenes = []
for i in oncoSites: # provides total oncogene region (may underestimate total size)
    i = i / 1000
    oncoGenes.append(i)
TIIIRegions = ['115227', '229041', '110541', '112997', '121167', '123547', '124428', '1397', '2126', '2390', '2593', '11486', '92527', '73379', '82455', '85949']

outFile = outputDir + '/chanceOncMut.txt'
output = open(outFile, 'w')
#write headers
output.write('Sample\tOncMutProb\tOncGeneProb\tTIIIProb\n')

perGeneOrSample = '1'
#Computes probability of any of the regional mutations occurring
if perGeneOrSample == '1':
    for x in samples:
        filename = inputDir + '/' + x + '/' + 'AF0_filtered.vcf'
        with open(filename, 'r') as target:
            oncoTotal = 0
            nonOncTotal = 0
            TIIITotal = 0
            otherTotal = 0
            for line in target:
                if '#' not in line and 'chr' in line: #skip the damn info
                    loc = line.split()[1]

                    #for oncogene mutations
                    #ask if oncogenic mutation and if not ask if in oncogene
                    if int(loc) in oncoSites: #is the var oncogenic?
                        oncoTotal += 1
                    elif int(loc)/1000 in oncoGenes: #is the var nonOncogenic exomic?
                        nonOncTotal += 1
                    elif str(int(loc)/1000) in TIIIRegions: #is the var TIII?
                        TIIITotal += 1

        #32 oncogenic sites
        oncoProb = float(oncoTotal)/32
        #15 probes covering 150bp - 32 bp that are oncogenic
        oncoGenProb = float(nonOncTotal)/(15 * 150 - 32)
        #17 probes covering 150bp
        TIIIProb = float(TIIITotal)/(17 * 150)

        #write to file
        output.write('%s\t%f\t%f\t%f\n' % (x, oncoProb, oncoGenProb, TIIIProb))

    target.close()

#The idea here is to look at individual oncogenes/regions across multiple samples and ask
#whether or not oncogenes are more preferentially mutated across samples
#completely unused right now
elif perGeneOrSample == '2':
    oncogeneDict = {'JAK2':[5073770], 'P53':[7577539, 7578402, 7577119], 'NRAS':[115256529, 115258747, 115258744], 'HRAS':[534287, 534288, 534289], 'KRAS':[25398284, 25380275], 'TET2':[106197266, 106197267, 106197268, 106197269, 106155172, 106155173, 106155174, 106155175], 'DNMT3A':[25457242, 25457243], 'IDH1':[209113112, 209113113], 'IDH2':[90631934, 90631838], 'GATA1':[48649700]}
