#The purpose of this script is to understand how well mutations in the
#dilution samples are found from the donor DNA.
#This script looks at variants that are AF=1 of the donor-in
#sample unfortunately there is no way to know if it unique because the seq
#of each sample alone failed.
import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--indir', '-i', required=True, type=str, help='Specifies the input directory containing all folders containing output analysis from a fermi analysis run.')
parser.add_argument('--donor', '-d', required=True, type=str, help='Name of the directory containing fermi analysis of donor dilution sample')
parser.add_argument('--recipients', '-r', required=True, type=str, help='Name of the directories containing fermi analysis of recipient dilution samples')


inputDir = args.indir
donorSample = args.donor
samples = args.recipients # this is a list
donorFile = 'AF1_filtered.vcf'
recipientFile = 'AF0_filtered.vcf'
recipientCheck = 'AF1_filtered.vcf' # to check if donor var is in recipient

consName = '_finalOutput.vcf'

spikeList = []
List24 = []
List25 = []
List28 = []
List29 = []
listOfLists = [List24, List25, List28, List29]

target = open(inputDir + donorSample + consName, 'r')

for line in target:
    if '#' not in line and 'chr' in line: #skip the damn info
        AO = line.split(';')[5]
        DP = line.split(';')[7]
        AF = line.split(';')[3]

        AONum = int(AO.split(',')[0][3:])
        AFNum = float(AF.split(',')[0][3:])
        DPNum = int(DP.split(',')[0][3:])

#get spike-in vars
        if AONum > 5 and (AFNum == 1 or AFNum == 0.5) and DPNum > 50:
            loc = line.split('\t')[1]
            spikeList.append(loc)

#get all variants
for i in listOfLists:
    index = listOfLists.index(i)
    target = open(inputDir + samples[index] + consName, 'r')
    for line in target:
        if '#' not in line and 'chr' in line: #skip the damn info
            loc = line.split('\t')[1]
            i.append(loc)

outFile = open(inputDir + 'DilutionEfficiency.txt', 'w')
numVars = len(spikeList)

for sample in listOfLists:
    counter = 0
    index = listOfLists.index(sample)
    for loc in sample:
        if loc in spikeList:
            counter += 1
    fraction = float(counter) / numVars * 100
    outFile.write('In sample %s %f percent of donor in variants were observed.\n' % (samples[index], fraction))

print len(spikeList)


