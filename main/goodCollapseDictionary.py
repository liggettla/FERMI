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
            #umi_seq = line[0:11]+line[-7:] #Abs dist from start/end compatible with miSeq/hiSeq
            umi_seq = line[0:11]+line.rstrip('\n')[-11:] #Abs dist from start/end compatible with miSeq/hiSeq
            umi_seq = umi_seq.rstrip('\n')
            read_seq = line.rstrip('\n')[6:-6]
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

            # it is important for duplex collapsing to make sure reads are of the same length
            # when not duplex collapsing this should always be true
            if not is_unique and len(sequences[umi_seq][0][0]) == len(read_seq):
                sequences[umi_seq][0].append(read_seq)
            elif is_unique:
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
    desiredReadLength = readLength # length specified by the user ie 120bp

    for umi in sequences:
        isReadGood = True
        finalRead = ''
        numReads = len(sequences[umi][0])
        header = sequences[umi][1]
        header = ''.join(header)
        quality = sequences[umi][2]
        quality = ''.join(quality)

        # to handle duplex collapse make sure overlap is larger than specified readlength
        # otherise just use max overlap region
        # if larger than desired length, then used desired length
        readLength = len(sequences[umi][0][0])
        if readLength > desiredReadLength:
            readLength = desiredReadLength

        # This is for troubleshooting, delete when done
        from pprint import pprint
        pprint(sequences[umi][0])

        # check that all sequences are of same length
        goodSet = checkLengths(sequences[umi][0])

        if goodSet:
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

######################################
# Check all captures are same length #
######################################
# maybe it would be better to take the most common length
# for ease just using first read at the moment
def checkLengths(umis):
    length = len(umis[0])
    goodSet = True
    for i in umis:
        if len(i) != length:
            goodSet = False

    return goodSet

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

    ########################
    # Collapse Unique UMIs #
    ########################
    # bin unique umis and derive consensus read based on input parameters of error tolerance
    for umi in sequences:
        # parse through the sequences defaultdict
        isReadGood = True
        finalRead = ''
        numReads = len(sequences[umi][0])
        header = sequences[umi][1]
        header = ''.join(header)
        quality = sequences[umi][2]
        quality = ''.join(quality)
        trimmed_quality = quality[0:readLength] # make quality string match read length

        # iterate through loci compute consensus base
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

                # calculate rate at which bases differ from consensus
                if errorRate == 'Y' and calcError:
                    rate = calcErrorRate(A,T,G,C)
                    errorRateList.append(rate)

        if isReadGood and numReads >= supportingReads:
            # if read checks out, record its coverage
            coverageList.append(covCounter)
            # add info to dataframe
            duplexDict[umi] = {'header':header, 'seq':finalRead, 'qual':trimmed_quality}
    return coverageList, errorRateList, duplexDict

def get_one_bp_mismatches(umi):
    mismatches = []
    bases = ['A','T','G','C']
    for e,i in enumerate(umi):
        for b in bases:
            mismatches.append(umi[:e] + b + umi[e+1:])
    mismatches=set(mismatches) # only includes unique strings
    return mismatches

def find_complementary_umis(duplexDict):
    from Bio.Seq import Seq
    # deDuplexList contains the sorted valid pairs for collapsing
    # in the form: deDuplexList[dictPairList[{},{}]]
    deDuplexList = []
    pairList = []

    #from pprint import pprint
    #pprint(duplexDict)

    for umi in duplexDict:
        if umi not in pairList:
            dictPairList = []
            tempList = []
            tempList.append(umi)
            oneMismatchList = get_one_bp_mismatches(umi) # get all possible mismatches for this umi
            for i in oneMismatchList:
                complement = str(Seq(i).reverse_complement()) # just search duplexDict for possible mismatches
                if complement in duplexDict:
                    tempList.append(complement)
            if len(tempList) == 2: # if only a pair then match then add to final dictionary
                pairList.extend(tempList) # avoid repeat scanning
                # record pair for collapsing
                dictPairList.append(duplexDict[tempList[0]])
                dictPairList.append(duplexDict[tempList[1]])
                deDuplexList.append(dictPairList)

    return deDuplexList

def collapse_paired_reads(deDuplexList, readLength, final_output_file):
    from Bio.Seq import Seq
    for pairedList in deDuplexList:
        seq1 = pairedList[0]['seq']
        seq2 = str(Seq(pairedList[1]['seq']).reverse_complement())
        finalRead = ''.join(
                seq1[base] if seq1[base] == seq2[base] else 'N'
                for base in range(readLength)
                )

        with open(final_output_file, 'a') as target:
            target.write('%s\n%s\n+\n%s\n' % (pairedList[0]['header'], finalRead, pairedList[0]['qual']))

def calcCoverageError(coverageList, errorRateList, errorRate):
    averageCoverage = mean(coverageList)

    if errorRate == 'Y': # clunky but passes to outputCov
        averageErrorRate = mean(errorRateList)
        return averageErrorRate, averageCoverage
    else:
        return 0, averageCoverage

def main():
    pass
