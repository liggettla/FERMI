import collections
from pprint import pprint

target = open('16_finalOutput.vcf', 'r')
output = open('AONum.txt', 'w')
AOList = []
infoDict = {}

count = 1

for line in target:
    if count > 54:
        AO = line.split(';')[5]
        DP = line.split(';')[7]
        AF = line.split(';')[3]
#these just trim off the AO= or DP= or AF= part of the string
        AONum = int(AO.split(',')[0][3:])
        AFNum = float(AF.split(',')[0][3:])
        DPNum = str(DP.split(',')[0][3:])
        if int(DPNum) > 50:
            output.write(DPNum + '\n')


    count+=1

orderedInfoDict = collections.OrderedDict(sorted(infoDict.items()))
#pprint(orderedInfoDict, output)
