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
            umi_seq = line[0:6]+line[-6:] #Abs dist from start/end compatible with miSeq/hiSeq
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

'''
This is an improvement on the original dictionary containing a lamda defined group of
lists.
This method should provide the same result as the previous lamdba dictionary but
should do so more quickly and in a more easily understood manner.
'''

def buildNestedDict(input_file, distance_stringency, pickleOut):
    #Dict format:
    #sortedSeqs = {'UMI':{'header': header, 'quality': quality, 'seqs':[seq1, seq2]}}
    sortedSeqs = {}
    target = open(input_file, 'r')
    umi_list = []
    position = 1
    is_unique = True

    #import pdb
    #pdb.set_trace()
    for line in target:
        if position == 1:
            header = line.rstrip('\n')
            position += 1
        elif position == 2:
            #Assumes UMI is flanking first and last 6bp of read
            umi_seq = line[0:6]+line[-6:] #Abs dist from start/end
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
                sortedSeqs[umi_seq] = {'header': header, 'quality': quality, 'seqs': [read_seq]}
            else:
                is_unique = True
                for umi in umi_list:
                    if is_unique:
                        if distance(umi_seq, umi) <= distance_stringency:
                            is_unique = False
                            umi_seq = umi

                #add to sortedSeqs database
                if is_unique:
                    umi_list.append(umi_seq)
                    sortedSeqs[umi_seq] = {'header': header, 'quality': quality, 'seqs': [read_seq]}
                elif not is_unique:
                    sortedSeqs[umi_seq]['seqs'].append(read_seq)
    target.close()

    ###########################
    #Write Full Data Structure#
    ###########################
    '''
    pickleFile = open(pickleOut, 'wb')
    pickle.dump(sortedSeqs, pickleFile)
    pickleFile.close()
    '''

    #import pprint
    #pprint.pprint(sortedSeqs)
    return sortedSeqs

#Collapses reads that come as a nested list dictionary
def collapseReadsListDict(sequences, varThresh, final_output_file, supportingReads, readLength):
    plus = '+'
    target = open(final_output_file, 'w')

    for umi in sequences:
        isReadGood = True
        finalRead = ''
        numReads = len(sequences[umi][0])
        header = sequences[umi][1]
        header = ''.join(header)
        quality = sequences[umi][2]
        quality = ''.join(quality)

        for base in range(readLength):
            if isReadGood:
                A = 0
                T = 0
                G = 0
                C = 0

                for seq in sequences[umi][0]:
                    if seq[base] == 'A':
                        A += 1
                    elif seq[base] == 'T':
                        T += 1
                    elif seq[base] == 'G':
                        G += 1
                    elif seq[base] == 'C':
                        C += 1

                if float(A)/numReads >= varThresh:
                    finalRead += 'A'
                elif float(T)/numReads >= varThresh:
                    finalRead += 'T'
                elif float(G)/numReads >= varThresh:
                    finalRead += 'G'
                elif float(C)/numReads >= varThresh:
                    finalRead += 'C'
                else:
                    #if too many errors ignore reads
                    isReadGood = False

        if isReadGood and numReads >= supportingReads:
            target = open(final_output_file, 'a')
            trimmed_quality = quality[6:-6] #MiSeq run
            target.write(header + '\n' + finalRead + '\n' + plus + '\n' + trimmed_quality + '\n')

#Collapses reads that are sorted by buildNestedDict()
def collapseNestedDict(sequences, varThresh, final_output_file, supportingReads, readLength):

    plus = '+'
    target = open(final_output_file, 'w')

    for umi in sequences:
        isReadGood = True
        finalRead = ''
        numReads = len(sequences[umi]['seqs'])
        header = sequences[umi]['header']
        quality = sequences[umi]['quality']

        for base in range(readLength):
            if isReadGood:
                A = 0
                T = 0
                G = 0
                C = 0

                for seq in sequences[umi]['seqs']:
                    if seq[base] == 'A':
                        A += 1
                    elif seq[base] == 'T':
                        T += 1
                    elif seq[base] == 'G':
                        G += 1
                    elif seq[base] == 'C':
                        C += 1

                if float(A)/numReads >= varThresh:
                    finalRead += 'A'
                elif float(T)/numReads >= varThresh:
                    finalRead += 'T'
                elif float(G)/numReads >= varThresh:
                    finalRead += 'G'
                elif float(C)/numReads >= varThresh:
                    finalRead += 'C'
                else:
                    #if too many errors ignore reads
                    isReadGood = False

        if isReadGood and numReads >= supportingReads:
            target = open(final_output_file, 'a')
            trimmed_quality = quality[6:-6]
            target.write(header + '\n' + finalRead + '\n' + plus + '\n' + trimmed_quality + '\n')

#Outdated and should be unused
def outputCov(distance_stringency):
    target = open(coverage_file, 'w')

    inputLines = sum(1 for line in open(input_file))
    outputLines = sum(1 for line in open(final_output_file))
    target.write("Eliminating only exact and close UMI matches:\n")
    target.write("# of Allowed UMI Mismatches: %d\n" %(distance_stringency))
    if not inputLines/4 == 0:
        target.write("Total # of Original UMIs: %d\n" % (inputLines/4))
        target.write("# of Unique UMIs: %d\n" % (outputLines/4))
    #avgCov = float(inputLines)/outputLines
    #target.write("Avg UMI Coverage: %r\n" % (avgCov))

    target.close()

def main():
    outputCov(distance_stringency)
