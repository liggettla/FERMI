#The purpose of this script is to take the output from oncogeneCoverage.sh
#and calculate a total coverage of the probed oncogenes

inputDir = '/media/alex/Extra/Dropbox/Code/fermiData/vcfs/'
inputFiles=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
outputFile = open('totalCoverageOut.txt', 'w')

for i in inputFiles:
    i = inputDir + str(i) + '_cov.txt'
    target = open(i, 'r')
    coverageList = []
    totalCoverage = 0

    for line in target:
        if line not in coverageList:
            totalCoverage += int(line)
            coverageList.append(totalCoverage)

    outputFile.write('Coverage for ' + i + ' is: ' + str(totalCoverage) + '\n')
    target.close()
outputFile.close()




