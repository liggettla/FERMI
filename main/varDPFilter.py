#################
#Filter Variants#
#################
# The purpose of this script is to filter the variants output in the final
# vcf files after running FERMI or the original TrueSeq analysis
# filters final vcf file to output either only AF=0 reads
# or AF=0.5 and AF=1 reads
import collections
from pprint import pprint

def vcfFilter(inDir, outputDir, inputFile, AOCutoff, DPCuttoff):

    inFiles = [inputFile]

    for file in inFiles:
        vcfOutputAF0 = outputDir + '/AF0_filtered.vcf' # includes AF=0
        vcfOutputAF1 = outputDir + '/AF1_filtered.vcf' # includes AF=1,0.5
        totalFiltered = outputDir + '/total_filtered.vcf' # includes all good vars
        plottable0 = outputDir + '/AF0_plottable.txt'
        plottable1 = outputDir + '/AF1_plottable.txt'
        totalPlottable = outputDir + '/total_plottable.txt'

        target = open(inputFile, 'r')
        vcfOut0 = open(vcfOutputAF0, 'w')
        vcfOut1 = open(vcfOutputAF1, 'w')
        vcfTotal = open(totalFiltered, 'w')
        dataOut0 = open(plottable0, 'w')
        dataOut1 = open(plottable1, 'w')
        dataTotal = open(totalPlottable, 'w')

        dataOut0.write('AO\tDP\tAF\n')
        dataOut1.write('AO\tDP\tAF\n')
        dataTotal.write('AO\tDP\tAF\n')


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
                if int(DPNum) > DPCuttoff:
                    if AONum > AOCutoff:
                        if AFNum == 0:
                            vcfOut0.write(line)
                            dataOut0.write(str(AONum) + '\t' + str(DPNum) + '\t' + str(AFNum) + '\n')
                            vcfTotal.write(line)
                            dataTotal.write(str(AONum) + '\t' + str(DPNum) + '\t' + str(AFNum) + '\n')
                        elif AFNum == 1 or AFNum == 0.5:
                            vcfOut1.write(line)
                            dataOut1.write(str(AONum) + '\t' + str(DPNum) + '\t' + str(AFNum) + '\n')
                            vcfTotal.write(line)
                            dataTotal.write(str(AONum) + '\t' + str(DPNum) + '\t' + str(AFNum) + '\n')

            else: #write the info lines
                vcfOut0.write(line)
                vcfOut1.write(line)
                vcfTotal.write(line)

        target.close()
        vcfOut0.close()
        vcfOut1.close()
        vcfTotal.close()
        dataOut0.close()
        dataOut1.close()
        dataTotal.close()

        #orderedInfoDict = collections.OrderedDict(sorted(infoDict.items()))
#pprint(orderedInfoDict, output)
