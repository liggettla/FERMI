#The purpose of the is script is to test the possibility
#that collapsing before aligning is causing problems (ITS NOT)

inDir = './testInput/'
outDir = './testOutput/'
r1Name = 'R1'
r2Name = 'R2'

r1In = open(inDir + r1Name + '.fastq', 'r')
r2In = open(inDir + r2Name + '.fastq', 'r')
r1Out = open(outDir + r1Name + '.fastq', 'w')
r2Out = open(outDir + r2Name + '.fastq', 'w')

count = 1
plus = '+'
for line in r1In:
    if count == 1:
        header = line.rstrip('\n')
        count += 1
    elif count == 2:
        trimmedSeq = line[6:-6]
        umiSeq = line[0:6]+line[-6:]
        umiSeq = umiSeq.rstrip('\n')
        count += 1
    elif count == 3:
        count += 1
    elif count == 4:
        trimmedQual = line[6:-6]
        count += 0
        r1Out.write(header + ':' + umiSeq + '\n' + trimmedSeq + '\n' + plus + '\n' + trimmedQual + '\n')


