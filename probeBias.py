'''
The purpose of this script is to look at probe bias and uses the data
coming from probeBias.sh.
'''

inputDir = '/media/alex/Extra/Dropbox/Code/fermiData/probeBias'
outputFile = '/media/alex/Extra/Dropbox/Code/fermiData/probeBias/AverageCoverage.txt'
output = open(outputFile, 'w')

#These are the oncogene probe sites
files=('50737', '75775', '75784', '75771', '1152565', '1152587', '1152587', '53428', '53428', '53428', '253982', '253802', '254572', '254572', '2091131', '2091131', '906319', '906318', '486497', '1982668')

for probe in files:
    inputFile = inputDir + '/' + probe + '.txt'
    target = open(inputFile, 'r')

    lastLine = 0
    total = 0
    count = 0
    for line in target:
        if line != lastLine:
            total += line
            count += 1
            lastLine = line

    average = total / count
    target.close()
    output.write('For Probe ' + probe + ' Coverage = ' + average)

output.close()












