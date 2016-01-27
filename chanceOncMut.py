#This purpose of this script is to ask the question of whether or not
#oncogenic mutations are being seen more often chance would predict based on
#the frequency of mutations in TIII and non-oncogenic regions as compared
#with the frequency seen at oncogenic locations
oncoSites = [5073770, 7577539,7577119,115256529,115258747,115258744,534287,534288,534289,25398284,25380275,106197266,106197267,106197268,106197269,106155172,106155173,106155174,25457242,25457243,209113112,209113113,90631934,90631838,48649700]
oncoGenes = []
for i in oncoSites:
    i = i / 1000
    oncoGenes.append(i)
TIIIRegions = ['115227', '229041', '110541', '112997', '121167', '123547', '124428', '1397', '2126', '2390', '2593', '11486', '92527', '73379', '82455', '85949']

output = open('../fermiData/chanceOncMut.txt', 'w')
#write headers
output.write('Sample\tOncMutProb\tOncGeneProb\tTIIIProb\n')

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
                nonOncTotal = 0
                TIIITotal = 0
                otherTotal = 0
                for line in target:
                    #The annoyance of dealing with the vcf info lines
                    if count > 54:
                        loc = line.split()[1]

                        #for oncogene mutations
                        #ask if oncogenic mutation and if not ask if in oncogene
                        if int(loc) in oncoSites: #is the var oncogenic?
                            oncoTotal += 1
                        elif int(loc)/1000 in oncoGenes: #is the var nonOncogenic exomic?
                            nonOncTotal += 1
                        elif str(int(loc)/1000) in TIIIRegions: #is the var TIII?
                            TIIITotal += 1

                    count += 1

            #32 oncogenic sites
            oncoProb = float(oncoTotal)/32
            #26 probes covering 150bp - 32 bp that are oncogenic
            oncoGenProb = float(nonOncTotal)/(26 * 150 - 32)
            #17 probes covering 150bp
            TIIIProb = float(TIIITotal)/(17 * 150)

            #write to file
            output.write('%i\t%f\t%f\t%f\n' % (filecount, oncoProb, oncoGenProb, TIIIProb))
            filecount += 1

        #skip these filenames b/c they don't exist
        elif filecount == 3 or filecount == 4:
            filecount += 1
        target.close()

#The idea here is to look at individual oncogenes/regions across multiple samples and ask
#whether or not oncogenes are more preferentially mutated across samples
#completely unused right now
elif perGeneOrSample == '2':
    oncogeneDict = {'JAK2':[5073770], 'P53':[7577539, 7578402, 7577119], 'NRAS':[115256529, 115258747, 115258744], 'HRAS':[534287, 534288, 534289], 'KRAS':[25398284, 25380275], 'TET2':[106197266, 106197267, 106197268, 106197269, 106155172, 106155173, 106155174, 106155175], 'DNMT3A':[25457242, 25457243], 'IDH1':[209113112, 209113113], 'IDH2':[90631934, 90631838], 'GATA1':[48649700]}
