#The purpose of this script is to compare the dilution samples in the
#1.5.2015 sequencing analysis or similar dilutions to understand how often
#a given mutation shows up when capturing and sequencing the same source
#DNA multiple times

count = 1
with open('18_finalOutput.vcf', 'r') as target:
    for line in target:
        #The annoyance of dealing with the vcf info lines
        if count > 54:
            print line.split()[1]
        count += 1
target.close()


