#The purpose of this script is to filter the variants output in the final vcf files
#after running FERMI or the original TrueSeq analysis
import collections
from pprint import pprint

inDir = './dataFiles/'
#inFiles = ['10_finalOutput', '13_finalOutput', '18_finalOutput', '20_finalOutput', '6_finalOutput', '9_finalOutput', '11_finalOutput', '14_finalOutput', '16_finalOutput', '19_finalOutput', '2_finalOutput', '7_finalOutput', '12_finalOutput', '15_finalOutput', '17_finalOutput', '1_finalOutput', '5_finalOutput', '8_finalOutput']
inFiles = ['21_finalOutput', '24_finalOutput', '25_finalOutput', '26_finalOutput', '27_finalOutput', '28_finalOutput', '29_finalOutput', '30_finalOutput']

for file in inFiles:
    input = inDir + file + '.vcf'
    vcfOutput = inDir + file + '_filtered.vcf'
    dataOutput = inDir + file + '_plottableData.txt'

    target = open(input, 'r')
    vcfOut = open(vcfOutput, 'w')
    dataOut = open(dataOutput, 'w')

    dataOut.write('AO\tDP\tAF\n')

    count = 1
    for line in target:
        if '#' not in line and 'chr' in line: #skip the damn info
            AO = line.split(';')[5]
            DP = line.split(';')[7]
            AF = line.split(';')[3]


#these just trim off the AO= or DP= or AF= part of the string
            AONum = int(AO.split(',')[0][3:])
            AFNum = float(AF.split(',')[0][3:])
            DPNum = str(DP.split(',')[0][3:])

#filtering parameters
            if int(DPNum) > 500:
                if AFNum == 0:
                    if AONum > 5:
                        vcfOut.write(line)
                        dataOut.write(str(AONum) + '\t' + str(DPNum) + '\t' + str(AFNum) + '\n')


        else: #write the info lines
            vcfOut.write(line)

    target.close()
    vcfOut.close()
    dataOut.close()

    #orderedInfoDict = collections.OrderedDict(sorted(infoDict.items()))
#pprint(orderedInfoDict, output)


'''
        if count <= 54:
            vcfOut.write(line)
        elif count > 54:
            AO = line.split(';')[5]
            DP = line.split(';')[7]
            AF = line.split(';')[3]
'''
