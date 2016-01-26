#This purpose of this script is to ask the question of whether or not
#oncogenic mutations are being seen more often chance would predict based on
#the frequency of mutations in TIII and non-oncogenic regions as compared
#with the frequency seen at oncogenic locations
oncoRegions = [5073770, 7577539,7577119,115256529,115258747,115258744,534287,534288,534289,25398284,25380275,106197266,106197267,106197268,106197269,106155172,106155173,106155174,25457242,25457243,209113112,209113113,90631934,90631838,48649700]
TIIIRegions = ['115227', '229041', '110541', '112997', '121167', '123547', '124428', '1397', '2126', '2390', '2593', '11486', '92527', '73379', '82455', '85949']

output = open('../fermiData/chanceOncMut.txt', 'w')

perGeneOrSample = raw_input('Analyse per sample or per gene (1/2): ')

#This needs to be improved, because right now I am manually switching between
#accurate TIII computation and Oncogene computation b/c not enough time to make automatic
if perGeneOrSample == '1':
    filecount = 1
    while filecount < 21:
        if filecount != 3 and filecount != 4: #no files 3 or 4 exist
            filename = '../fermiData/vcfs/' + str(filecount) + '_finalOutput.vcf'
            count = 1
            with open(filename, 'r') as target:
                oncoTotal = 0
                otherTotal = 0
                for line in target:
                    #The annoyance of dealing with the vcf info lines
                    if count > 54:
                        loc = line.split()[1]

                        #for oncogenic mutations
                        '''
                        if int(loc) in oncoRegions: #is the var oncogenic?
                            oncoTotal += 1
                        else: #is the var nonOncogenic?
                            otherTotal += 1
                        '''

                        #for TIII mutations
                        #This cuts off the last 3 digits of each location to match the trimmed
                        #list locations
                        if str(int(loc)/1000) in TIIIRegions: #is the var oncogenic?
                            otherTotal += 1
                        else: #is the var nonOncogenic?
                            oncoTotal += 1

                    count += 1

            #26 oncogenic sites
            oncoProb = float(oncoTotal)/26
            #43 probes covering 150bp - 26 bp that are oncogenic
            otherProb = float(otherTotal)/(43 * 150 - 26)

            #different ways to write the data to file
            #output.write('%i\t%f\t%f\n' % (filecount, oncoProb, otherProb))
            output.write('Sample %i \nOnc Mut Prob: %f \nNon-Onc Mut Prob: %f\n' % (filecount, oncoProb, otherProb))
            filecount += 1

        #skip these filenames b/c they don't exist
        elif filecount == 3 or filecount == 4:
            filecount += 1
        target.close()

#The idea here is to look at individual oncogenes/regions across multiple samples and ask
#whether or not oncogenes are more preferentially mutated across samples
elif perGeneOrSample == '2':
    oncogeneDict = {'JAK2':[5073770], 'P53':[7577539, 7578402, 7577119], 'NRAS':[115256529, 115258747, 115258744], 'HRAS':[534287, 534288, 534289], 'KRAS':[25398284, 25380275], 'TET2':[106197266, 106197267, 106197268, 106197269, 106155172, 106155173, 106155174, 106155175], 'DNMT3A':[25457242, 25457243], 'IDH1':[209113112, 209113113], 'IDH2':[90631934, 90631838], 'GATA1':[48649700]}
#Iterate through each of the genes
    for gene in oncogeneDict:

#Iterate through each of the files
        while filecount < 21:
            if filecount != 3 and filecount != 4: #no files 3 or 4 exist
                filename = '../fermiData/vcfs/' + str(filecount) + '_finalOutput.vcf'
                count = 1
                with open(filename, 'r') as target:






