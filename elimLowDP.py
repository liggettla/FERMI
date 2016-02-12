#This purpose of this script is to eliminate those variants that are found at low DP
#This also eliminates germline variants

inputFile = open('16_finalOutput.vcf', 'r')
outputFile = open('16_finalOutput_highDP.vcf', 'w')

count = 1
for line in inputFile:
    if count > 54:
        DP = line.split(';')[7]
        AF = line.split(';')[3]
        AO = line.split(';')[5]
        DPNum = str(DP.split(',')[0][3:])
        AFNum = float(AF.split(',')[0][3:])
        AONum = int(AO.split(',')[0][3:])

        if int(DPNum) > 500:
            if AFNum == 0:
                if AONum > 5:
                    outputFile.write(line)
    else:
        outputFile.write(line)

    count += 1
