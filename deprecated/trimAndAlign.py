#The purpose of this script is to generate fastq file without UMI sequences
#and then align and variant call. This can be useful for finding true VAF in donor
#samples for FERMI testing for example.

#files = ['21_ATCACGAC_L008_R1_001.fastq', '25_ACCCAGCA_L008_R1_001.fastq', '27_CCCAACCT_L008_R1_001.fastq', '29_GAAACCCA_L008_R1_001.fastq', '24_ACAAACGG_L008_R1_001.fastq', '26_AACCCCTC_L008_R1_001.fastq', '28_CACCACAC_L008_R1_001.fastq', '30_TGTGACCA_L008_R1_001.fastq']
from os import system

files = ['21_ATCACGAC_L008_R1_001.fastq']

for i in files:
    target = open(i, 'r')
    outName = i.strip('.fastq') + '_noUMI.fastq'
    outFile = open(outName, 'w')

    plus = '+'
    position = 1
    for line in target:
        if position == 1:
            header = line.rstrip('\n')
            position += 1
        elif position == 2:
            read_seq = line[6:] #trims off 5' UMI
            position += 1
        elif position == 3:
            position += 1
        elif position == 4:
            quality = line[6:]
            position = 1
            outFile.write(header + '\n' + read_seq + plus + '\n' + quality)

    REF = '/vol3/home/liggettl/refgenomes/hg19.fa'

    bamOut = outName.strip('fastq') + 'bam'
    vcfOut = outName.strip('fastq') + 'vcf'

    system("bwa mem %s %s | samtools view -bS - | samtools sort > %s" % (REF, outName, bamOut))
    system("samtools index %s" % (bamOut))

    #Cluster
    system("/vol3/home/liggettl/TruSeqPanel/Scripts/freebayes/freebayes --fasta-reference %s %s > %s" % (REF, bamOut, vcfOut))
