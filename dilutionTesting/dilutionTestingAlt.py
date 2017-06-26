#!/usr/bin/env python

#The purpose of this script is to understand how well mutations in the
#dilution samples are found from the donor DNA.
#This script looks at variants that are AF=1 of the donor-in
#sample unfortunately there is no way to know if it unique because the seq
#of each sample alone failed.
import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--indir', '-i', required=True, type=str, help='Specifies the input directory containing all folders containing output analysis from a fermi analysis run.')
parser.add_argument('--donor', '-d', required=True, type=str, help='Name of the directory containing fermi analysis of donor dilution sample')
parser.add_argument('--controlRecipient', '-c', required=True, type=str, help='Name of the directory containing the blank recipient sample into which no donor DNA has been input')
parser.add_argument('--recipients', '-r', required=True, nargs='*', type=str, help='Name of the directories containing fermi analysis of recipient dilution samples')
parser.add_argument('--nonunique', '-n', action='store_true', help='Instead of looking for variants that are unique to the donor sample, use variants that are homo/het in donor sample and more rare than 1/100 in the recipient samples.')

args = parser.parse_args()

#############
#Input Files#
#############
inputDir = args.indir
donorSample = args.donor # donor DNA
controlRec = args.controlRecipient # blank recipient DNA
recSamples = args.recipients # this is a list

# if nonunique rare vars from recipient files should still be included
if args.nonunique:
    nonUnique = 'Y'
elif not args.nonunique:
    nonUnique = 'n'
# default to using nonUniques
else:
    nonUnique = 'Y'


donorFile = inputDir + '/' + donorSample + '/AF1_filtered.vcf'
# to check if donor var is in recipient
recipientCheck = inputDir + '/' + donorSample + 'AF1_filtered.vcf'

#############################
# Build Accurate Donor List #
#############################
# data structure looks like this donorDict = {loc:{ref:base, alt:base}}
if nonUnique == 'Y':
    donorDict = {}
    target = open(donorFile, 'r')
    for line in target:
        if '#' not in line and 'chr' in line: #skip the damn info
            ref = line.split('\t')[3]
            alt = line.split('\t')[4]
            loc = line.split('\t')[1]
            AO = float(line.split(';')[5].strip('AO='))
            DP = float(line.split(';')[7].strip('DP='))

            donorDict[loc] = {'ref':ref, 'alt':alt, 'AO':AO, 'DP':DP}

########################
# Build Recipient List #
########################
if nonUnique == 'Y':
    recipientDict = {}
    control = inputDir + '/' + controlRec + '/onlyProbedRegions.vcf'
    target = open(control, 'r')

    for line in target:
        if '#' not in line and 'chr' in line: #skip the damn info
            ref = line.split('\t')[3]
            alt = line.split('\t')[4]
            loc = line.split('\t')[1]
            AO = float(line.split(';')[5].strip('AO='))
            DP = float(line.split(';')[7].strip('DP='))

            # only include vars that are more rare than 1:100 in final analysis
            if AO/DP > 0.01:
                recipientDict[loc] = {'ref':ref, 'alt':alt}

#########################
# Get Unique Donor Vars #
#########################
# get vars that are found at high freq in the donor sample and at low freq in the recipient sample
if nonUnique == 'Y':
    uniqDonor = {}
    for i in donorDict:
        if donorDict[i]['alt'] not in recipientDict:
            uniqDonor[i] = {'ref':donorDict[i]['ref'], 'alt':donorDict[i]['alt'], 'AO':donorDict[i]['AO'], 'DP':donorDict[i]['DP']}



# Find Vars in Recipients #
if nonUnique == 'Y':
    results = inputDir + '/dilutionResults.txt'
    resultsTarget = open(results, 'w')

    for sample in recSamples:
        sampleDict = {}
        x = inputDir + '/' + sample + '/total_filtered.vcf'
        target = open(x, 'r')

        for line in target:
            if '#' not in line and 'chr' in line: #skip the damn info
                ref = line.split('\t')[3]
                alt = line.split('\t')[4]
                loc = line.split('\t')[1]
                AO = float(line.split(';')[5].strip('AO='))
                DP = float(line.split(';')[7].strip('DP='))

                sampleDict[loc] = {'ref':ref, 'alt':alt, 'AO':AO, 'DP':DP}



                if loc in uniqDonor:
                    AF = AO/DP
                    resultsTarget.write('Sample: ' + str(sample) + '\n' + 'Location: ' + str(loc) + '\n' \
                            + 'AO: ' + str(sampleDict[loc]['AO']) + '\n' \
                            + 'DP: ' + str(sampleDict[loc]['DP']) + '\n' \
                            + 'AF: ' + str(AF) + '\n')
                    print sample
                    print loc,sampleDict[loc]





#################################
# Output AF in dilution samples #
#################################
if nonUnique == 'n':

    target.write('Sample: ' + str(i) + '\n' + 'Location: ' + str(loc) + '\n' \
            + 'AO: ' + str(AONum) + '\n' \
            + 'DP: ' + str(DPNum) + '\n' \
            + 'AF: ' + str(AFNum) + '\n')






######################
#Build Donor Loc List#
######################
donorList = []
target = open(donorFile, 'r')
for line in target:
    if '#' not in line and 'chr' in line: #skip the damn info
        loc = line.split('\t')[1]
        donorList.append(loc)
target.close()

##########################
#Build Recipient Loc Dict#
##########################
recAF0Dict = {}
recAF1Dict = {}
for i in recSamples:
    A1 = inputDir + '/' + i + '/AF1_filtered.vcf'
    target1 = open(A1, 'r')

    AF1List = []
    for line in target1:
        if '#' not in line and 'chr' in line: #skip the damn info
            loc = line.split('\t')[1]
            AF1List.append(loc)

    recAF1Dict[i] = AF1List

target1.close()

#####################
#Get Uniq Donor Vars#
#####################
commonDonor = []
uniqDonor = []
for i in donorList:
    for j in recAF1Dict:
        if i in recAF1Dict[j]:
            commonDonor.append(i)

for i in donorList:
    if not i in commonDonor:
        uniqDonor.append(i)

#####################
#Output Observations#
#####################
results = inputDir + '/dilutionResults.txt'
target = open(results, 'w')

for i in recSamples:
    A0 = inputDir + '/' + i + '/AF0_filtered.vcf'
    target0 = open(A0, 'r')

    for line in target0:
        if '#' not in line and 'chr' in line: #skip the damn info
            loc = line.split('\t')[1]
            AO = line.split(';')[5]
            DP = line.split(';')[7]
            #AF = line.split(';')[3]

            AONum = float(AO.split(',')[0][3:])
            DPNum = float(DP.split(',')[0][3:])
            AFNum = AONum / DPNum

            if loc in uniqDonor:
                target.write('Sample: ' + str(i) + '\n' + 'Location: ' + str(loc) + '\n' \
                        + 'AO: ' + str(AONum) + '\n' \
                        + 'DP: ' + str(DPNum) + '\n' \
                        + 'AF: ' + str(AFNum) + '\n')

target.close()
target0.close()

################
# Output Plots #
################
import numpy as np
import matplotlib.pyplot as plt



