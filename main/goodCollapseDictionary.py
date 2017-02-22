'''
Unless otherwise noted, this is the most up-to-date artificial consensus read
derivation script that should now be fully compatible with both MiSeq and HiSeq runs
and independent of amplicon length.

This is similar to identifyUniqueUMI.py and quickCollapse.py but attempts to
improve upon processing speed by first putting all reads into a dictionary
with UMIs as keys and sequences as sequences in corresponding lists

For a HiSeq run with 150 cycle chemistry to process correctly, concatenateUMI.py must
be run before the methods in this script.
'''

import pdb
from collections import defaultdict
from Levenshtein import distance #string length can be different
from itertools import islice
import pickle
from numpy import mean
from errorRate import calcErrorRate

####################
# Build Dictionary #
####################

#Tracks quality score and header as well as sequence and UMI
#by reading in info from input_file into a dict containing three nested lists
#associated with each UMI key
def buildListDict(input_file, distance_stringency, pickleOut):
    #Dict format: {'UMI_1' : (Seqs, First_Header, First_quality), 'UMI_2' : (Seqs, First_Header, First_quality)}
    sequences = defaultdict(lambda:([],[],[]))
    target = open(input_file, 'r')
    umi_list = []
    position = 1
    is_unique = True

    for line in target:
        if position == 1:
            header = line.rstrip('\n')
            position += 1
        elif position == 2:
            #Assumes UMI is flanking first and last 6bp of read
            umi_seq = line[0:11]+line[-6:] #Abs dist from start/end compatible with miSeq/hiSeq
            umi_seq = umi_seq.rstrip('\n')
            read_seq = line[6:-6]
            position += 1
        elif position == 3:
            position += 1
        elif position == 4:
            quality = line.rstrip('\n')
            position = 1

            if not bool(umi_list):
                umi_list.append(umi_seq)
            else:
                is_unique = True
                for umi in umi_list:
                    if is_unique:
                        if distance(umi_seq, umi) <= distance_stringency:
                            is_unique = False
                            umi_seq = umi

            sequences[umi_seq][0].append(read_seq)

            #check if header slot is empty prevents multiple entries error
            if is_unique and not bool(sequences[umi_seq][1]):
                sequences[umi_seq][1].append(header)
                sequences[umi_seq][2].append(quality)

    target.close()
    return sequences

#######################
# Collapse Dictionary #
#######################

#Collapses reads that come as a nested list dictionary
def collapseReadsListDict(sequences, varThresh, final_output_file, supportingReads, readLength, errorRate, badBaseSubstitute):
    errorRateList = []
    plus = '+'
    target = open(final_output_file, 'w')
    coverageList = [] # used to tally coverage for each UMI pair

    for umi in sequences:
        isReadGood = True
        finalRead = ''
        numReads = len(sequences[umi][0])
        header = sequences[umi][1]
        header = ''.join(header)
        quality = sequences[umi][2]
        quality = ''.join(quality)

        for base in range(readLength):
            covCounter = 0 # used for quality coverage calc

            if isReadGood:
                A = 0
                T = 0
                G = 0
                C = 0

                for seq in sequences[umi][0]:
                    covCounter += 1 # count number of supporting reads

                    if seq[base] == 'A':
                        A += 1
                    elif seq[base] == 'T':
                        T += 1
                    elif seq[base] == 'G':
                        G += 1
                    elif seq[base] == 'C':
                        C += 1

                calcError = True # used to include a base in error rate calc

                if float(A)/numReads >= varThresh:
                    finalRead += 'A'
                elif float(T)/numReads >= varThresh:
                    finalRead += 'T'
                elif float(G)/numReads >= varThresh:
                    finalRead += 'G'
                elif float(C)/numReads >= varThresh:
                    finalRead += 'C'
                else:
                    # if there are too many errors either substitute or invalidate
                    if badBaseSubstitute: # ignore base but keep read
                        isReadGood = True
                        finalRead += 'N'
                    else: # ignore entire read
                        isReadGood = False

                    calcError = False

                if errorRate == 'Y' and calcError:
                    rate = calcErrorRate(A,T,G,C)
                    errorRateList.append(rate)

        if isReadGood and numReads >= supportingReads:
            coverageList.append(covCounter) # if read checks out, record its coverage

            target = open(final_output_file, 'a')
            trimmed_quality = quality[6:-6]
            trimmed_quality = trimmed_quality[0:readLength] # make quality string match read length
            target.write(header + '\n' + finalRead + '\n' + plus + '\n' + trimmed_quality + '\n')

    averageCoverage = mean(coverageList)

    if errorRate == 'Y': # clunky but passes to outputCov
        averageErrorRate = mean(errorRateList)
        return averageErrorRate, averageCoverage
    else:
        return 0, averageCoverage

###################
# Duplex Collapse #
###################
def duplexCollapse(sequences, varThresh, final_output_file, supportingReads, readLength, errorRate, badBaseSubstitute):
    from Bio.Seq import Seq

    duplexDict = {}
    errorRateList = []
    plus = '+'
    target = open(final_output_file, 'w')
    coverageList = [] # used to tally coverage for each UMI pair

    for umi in sequences:
        isReadGood = True
        finalRead = ''
        numReads = len(sequences[umi][0])
        header = sequences[umi][1]
        header = ''.join(header)
        quality = sequences[umi][2]
        quality = ''.join(quality)
        trimmed_quality = quality[0:readLength] # make quality string match read length


        for base in range(readLength):
            covCounter = 0 # used for quality coverage calc

            if isReadGood:
                A = 0
                T = 0
                G = 0
                C = 0

                for seq in sequences[umi][0]:
                    covCounter += 1 # count number of supporting reads

                    if seq[base] == 'A':
                        A += 1
                    elif seq[base] == 'T':
                        T += 1
                    elif seq[base] == 'G':
                        G += 1
                    elif seq[base] == 'C':
                        C += 1

                calcError = True # used to include a base in error rate calc

                if float(A)/numReads >= varThresh:
                    finalRead += 'A'
                elif float(T)/numReads >= varThresh:
                    finalRead += 'T'
                elif float(G)/numReads >= varThresh:
                    finalRead += 'G'
                elif float(C)/numReads >= varThresh:
                    finalRead += 'C'
                else:
                    # if there are too many errors either substitute or invalidate
                    if badBaseSubstitute: # ignore base but keep read
                        isReadGood = True
                        finalRead += 'N'
                    else: # ignore entire read
                        isReadGood = False

                    calcError = False

                if errorRate == 'Y' and calcError:
                    rate = calcErrorRate(A,T,G,C)
                    errorRateList.append(rate)

        if isReadGood and numReads >= supportingReads:
            # if read checks out, record its coverage
            coverageList.append(covCounter)
            # add info to dataframe
            duplexDict[umi] = {'header':header, 'seq':finalRead, 'qual':trimmed_quality}

    ###########################
    # Find Complementary UMIs #
    ###########################
    deDuplexDict = {}
    for umi in duplexDict.keys():
        complement = str(Seq(umi).complement())
        if complement in duplexDict:
            deDuplexDict[umi] = duplexDict[umi]
            deDuplexDict[complement] = duplexDict[complement]
            del duplexDict[umi]
            del duplexDict[complement]
        else:
            tempDict = {}
            for i in duplexDict:
                if distance(complement, i) == 1:
                    tempDict[i]=duplexDict[i]

            if len(tempDict) == 1:
                deDuplexDict[umi] = duplexDict[umi]
                deDuplexDict[complement] = duplexDict[complement]
                del duplexDict[i]
                del duplexDict[umi]

    ###############################
    # Collapse Complementary UMIs #
    ###############################
    plus = '+'
    for umi in deDuplexDict.keys():
        finalRead = ''

        try:
            refRead = deDuplexDict[umi]['seq']
            complement = str(Seq(umi).complement())
            complementaryRead = str(Seq(deDuplexDict[complement]['seq']).complement())

            for base in range(readLength):
                if refRead[base] == complementaryRead[base]:
                    finalRead += refRead[base]
                else:
                    finalRead += 'N'

            target = open(final_output_file, 'a')
            target.write(deDuplexDict[umi]['header'] + '\n' + finalRead + '\n' + plus + '\n' + deDuplexDict[umi]['qual'] + '\n')

            del deDuplexDict[umi]
            del deDuplexDict[complement]

        except:
            pass

    averageCoverage = mean(coverageList)

    if errorRate == 'Y': # clunky but passes to outputCov
        averageErrorRate = mean(errorRateList)
        return averageErrorRate, averageCoverage
    else:
        return 0, averageCoverage

def main():
    pass
