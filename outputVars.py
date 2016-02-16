inputDir = './dataFiles/'

unfil8 = ['24_finalOutput.vcf', '25_finalOutput.vcf', '26_finalOutput.vcf', '27_finalOutput.vcf']
unfil9 = ['28_finalOutput.vcf', '29_finalOutput.vcf', '30_finalOutput.vcf']

fil8 = ['24_finalOutput_filtered.vcf', '25_finalOutput_filtered.vcf']
fil9 = ['28_finalOutput_filtered.vcf', '29_finalOutput_filtered.vcf']

dataSetList = [unfil8, unfil9, fil8, fil9]


out = inputDir + 'repeatability.txt'
outFile = open(out, 'w')

for group in dataSetList:
    for i in group:
        outFile.write('Starting with file %s \n' % (i))
        list1 = []

        file1 = open(inputDir + i, 'r')

        for line in file1:
            if '#' not in line and 'chr' in line:
                loc = line.split('\t')[1]
                list1.append(loc)

        for j in group:
            if i != j:
                list2 = []
                list3 = []
                outFile.write('Comparing to file %s \n' % (j))
                file2 = open(inputDir + j, 'r')
                for line in file2:
                    if '#' not in line and 'chr' in line:
                        loc = line.split('\t')[1]
                        list2.append(loc)
                for item in list1:
                    if item in list2:
                        list3.append(item)
                        outFile.write(item)
                        outFile.write('\n')
                fraction1 = str(float(len(list3))/len(list1) * 100)
                fraction2 = str(float(len(list3))/len(list2) * 100)
                outFile.write('Shared sequences are %s percent of file %s \n' % (fraction1, i))
                outFile.write('Shared sequences are %s percent of file %s \n' % (fraction2, j))

                file2.close()
        file1.close()




