#This purpose of this script is to ask the question of whether or not
#oncogenic mutations are being seen more often chance would predict based on
#the frequency of mutations in TIII and non-oncogenic regions as compared
#with the frequency seen at oncogenic locations
oncoList = [5073770, 7577539,7577119,115256529,115258747,115258744,534287,534288,534289,25398284,25380275,106197266,106197267,106197268,106197269,106155172,106155173,106155174,25457242,25457243,209113112,209113113,90631934,90631838,48649700]

oncList = []
nonOncList = []

filecount = 1
if filecount != 3 or filecount != 4: #no files 3 or 4 exist
    filename = '../fermiData/vcfs/' + filecount + '_finalOutput.vcf'
    filecount += 1
    count = 1
    with open(filename, 'r') as target:
        for line in target:
            #The annoyance of dealing with the vcf info lines
            if count > 54:
                i.append(line.split()[1])
            count += 1
    target.close()

#skip these filenames b/c they don't exist
elif filecount == 3 or filecount == 4:
    filecount += 1





