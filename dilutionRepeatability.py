#The purpose of this script is to compare the dilution samples in the
#1.5.2015 sequencing analysis or similar dilutions to understand how often
#a given mutation shows up when capturing and sequencing the same source
#DNA multiple times

list18 = []
list19 = []
list20 = []

with open('../fermiData/vcfs/18_finalOutput.vcf', 'r') as target:
    count = 1
    for line in target:
        #The annoyance of dealing with the vcf info lines
        if count > 54:
            if ';AF=0' not in line: #only look at homo/het vars
                list18.append(line.split()[1])
        count += 1
target.close()

with open('../fermiData/vcfs/19_finalOutput.vcf', 'r') as target:
    count = 1
    for line in target:
        #The annoyance of dealing with the vcf info lines
        if count > 54:
            if ';AF=0' not in line:
                list19.append(line.split()[1])
        count += 1
target.close()

with open('../fermiData/vcfs/20_finalOutput.vcf', 'r') as target:
    count = 1
    for line in target:
        #The annoyance of dealing with the vcf info lines
        if count > 54:
            if ';AF=0' not in line:
                list20.append(line.split()[1])
        count += 1
target.close()

#Compare lists
#18 and 19
recur = 0
uniq = 0
total = 0
for i in list18:
    if i in list19:
        recur += 1
    elif i not in list19:
        uniq += 1
total = recur + uniq
print ('Percent Common Vars between 18 and 19: ' + str(float(recur)/float(total)*100))

#18 and 20
recur = 0
uniq = 0
total = 0
for i in list18:
    if i in list20:
        recur += 1
    elif i not in list20:
        uniq += 1
total = recur + uniq
print ('Percent Common Vars between 18 and 20: ' + str(float(recur)/float(total)*100))

#19 and 20
recur = 0
uniq = 0
total = 0
for i in list19:
    if i in list20:
        recur += 1
    elif i not in list20:
        uniq += 1
total = recur + uniq
print ('Percent Common Vars between 19 and 20: ' + str(float(recur)/float(total)*100))
