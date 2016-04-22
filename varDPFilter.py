# The purpose of this script is to filter the variants output in the final
# vcf files after running FERMI or the original TrueSeq analysis
import collections
from pprint import pprint

def vcfFilter(inDir, outputDir, inputFile):

    inFiles = [inputFile]

    for file in inFiles:
        vcfOutputAF0 = outputDir + '/AF0_filtered.vcf' # includes AF=0
        vcfOutputAF1 = outputDir + '/AF1_filtered.vcf' # includes AF=1,0.5
        plottable0 = outputDir + '/AF0_plottable.txt'
        plottable1 = outputDir + '/AF1_plottable.txt'

        target = open(inputFile, 'r')
        vcfOut0 = open(vcfOutputAF0, 'w')
        vcfOut1 = open(vcfOutputAF1, 'w')
        dataOut0 = open(plottable0, 'w')
        dataOut1 = open(plottable1, 'w')

        dataOut0.write('AO\tDP\tAF\n')
        dataOut1.write('AO\tDP\tAF\n')

        count = 1
        for line in target:
            if '#' not in line and 'chr' in line: #skip the damn info
                AO = line.split(';')[5]
                DP = line.split(';')[7]
                AF = line.split(';')[3]

                # these just trim off the AO= or DP= or AF= part of the string
                AONum = int(AO.split(',')[0][3:])
                AFNum = float(AF.split(',')[0][3:])
                DPNum = str(DP.split(',')[0][3:])

                # filtering parameters
                if int(DPNum) > 500:
                    if AONum > 5:
                        if AFNum == 0:
                            vcfOut0.write(line)
                            dataOut0.write(str(AONum) + '\t' + str(DPNum) + '\t' + str(AFNum) + '\n')
                        elif AFNum == 1 or AFNum == 0.5:
                            vcfOut1.write(line)
                            dataOut1.write(str(AONum) + '\t' + str(DPNum) + '\t' + str(AFNum) + '\n')


            else: #write the info lines
                vcfOut0.write(line)
                vcfOut1.write(line)

        target.close()
        vcfOut0.close()
        vcfOut1.close()
        dataOut0.close()
        dataOut1.close()

        #orderedInfoDict = collections.OrderedDict(sorted(infoDict.items()))
#pprint(orderedInfoDict, output)
