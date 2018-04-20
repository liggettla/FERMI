#!/usr/bin/env bash

# laptop
dir="/media/alex/Extra"
# desktop
#dir="/media/alex/mainstorage/LinuxDropbox"

# 1b,c,d, 3a,c
#./outputcontrol -d $dir/Dropbox/Code/Sequencing/8.24.2017_OneSupportinReadPublicationFiles -r $dir/Dropbox/Code/ReferenceGenomes/hg19.fa -i 4r1.fastq 5r1.fastq 6r1.fastq 7r1.fastq 9r1.fastq 10r1.fastq 11r1.fastq 12r1.fastq 13r1.fastq 14r1.fastq 15r1.fastq 16r1.fastq 17r1.fastq 18r1.fastq -s allVariants.pkl -p allVariants.pkl

# 2c
#./plotRsquaredvsAge

# 2c?
#./oncogenicChanges -o $dir/Dropbox/Code/FERMI/repeatSequencingAnalysis -i $dir/Dropbox/Code/Sequencing/8.24.2017_OneSupportinReadPublicationFiles -s 1r1.fastq 3r1.fastq 4r1.fastq 5r1.fastq 6r1.fastq 7r1.fastq 9r1.fastq 10r1.fastq 11r1.fastq 12r1.fastq 13r1.fastq 14r1.fastq 15r1.fastq 16r1.fastq 17r1.fastq 18r1.fastq -k samplekey.txt

# 3b
#./quantifyCPG -i /media/alex/Extra/Dropbox/Code/Sequencing/8.24.2017_OneSupportinReadPublicationFiles/7r1.fastq/onlyProbedRegions.vcf -f 1 -d /media/alex/Extra/Dropbox/Code/Sequencing/8.24.2017_OneSupportinReadPublicationFiles/relabeledVCFs -m a.vcf b.vcf d.vcf e.vcf f.vcf g.vcf i.vcf j.vcf l.vcf m.vcf n.vcf o.vcf p.vcf q.vcf r.vcf s.vcf t.vcf u.vcf v.vcf -r /media/alex/Extra/Dropbox/Code/ReferenceGenomes/hg19.fa -pq

# 4a,b
#./multiSampleVAFRep.py -i /media/alex/Extra/Dropbox/Code/Sequencing/8.24.2017_OneSupportinReadPublicationFiles -o /media/alex/Extra/Dropbox/Code/FERMI/paperGeneration -s 15r1.fastq -p 15r1.fastq -r 0.3 -z

# 5b
./stoopidPlotting
