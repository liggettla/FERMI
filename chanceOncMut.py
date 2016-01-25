#This purpose of this script is to ask the question of whether or not
#oncogenic mutations are being seen more often chance would predict based on
#the frequency of mutations in TIII and non-oncogenic regions as compared
#with the frequency seen at oncogenic locations
oncoRegions = [5073770, 7577539,7577119,115256529,115258747,115258744,534287,534288,534289,25398284,25380275,106197266,106197267,106197268,106197269,106155172,106155173,106155174,25457242,25457243,209113112,209113113,90631934,90631838,48649700]

output = open('../fermiData/chanceOncMut.txt', 'w')

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

                    if loc in oncoRegions: #is the var oncogenic?
                        oncoTotal += 1
                    else: #is the var nonOncogenic?
                        otherTotal += 1
                count += 1

        #26 oncogenic sites
        oncoProb = float(oncoTotal)/26
        #43 probes covering 150bp - 26 bp that are oncogenic
        otherProb = float(otherTotal)/(43 * 150 - 26)
        output.write('Sample %i \nOnc Mut Prob: %f \nNon-Onc Mut Prob: %f\n' % (filecount, oncoProb, otherProb))
        filecount += 1

    #skip these filenames b/c they don't exist
    elif filecount == 3 or filecount == 4:
        filecount += 1





    target.close()
