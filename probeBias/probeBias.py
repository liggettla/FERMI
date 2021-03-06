from ..main import parseLine

'''
The purpose of this script is to look at probe bias and uses the data
coming from probeBias.sh.
'''

inputDir = '/media/alex/Extra/Dropbox/Code/fermiData/probeBias'
outputFile = '/media/alex/Extra/Dropbox/Code/fermiData/probeBias/AverageCoverage.txt'
output = open(outputFile, 'w')

#Oncogene probe sites
#files=('50737', '75775', '75784', '75771', '1152565', '1152587', '1152587', '53428', '53428', '53428', '253982', '253802', '254572', '254572', '2091131', '2091131', '906319', '906318', '486497', '1982668')
#TIII probe sites
files=('115227', '229041', '110541', '112997', '121167', '123547', '124428', '1397', '2126', '2390', '2593', '11486', '92527', '73379', '82455', '85949')

for probe in files:
    inputFile = inputDir + '/' + probe + '.txt'
    target = open(inputFile, 'r')

    lastLine = 0
    total = 0
    count = 0
    for line in target:
        line = int(line)
        if line != lastLine:
            total += line
            count += 1
            lastLine = line

    average = total / count
    target.close()
    output.write('For Probe ' + str(probe) + ' Coverage = ' + str(average) + '\n')

output.close()












