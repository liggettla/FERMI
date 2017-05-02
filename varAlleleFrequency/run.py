#!/usr/bin/env python

from os import system

sampleList = ('4r1.fastq', '5r1.fastq', '6r1.fastq', '7r1.fastq', '8r1.fastq', '9r1.fastq', '10r1.fastq', '11r1.fastq', '12r1.fastq', '13r1.fastq', '14r1.fastq', '15r1.fastq', '16r1.fastq', '17r1.fastq', '18r1.fastq', '19r1.fastq', '20r1.fastq', '21r1.fastq', '22r1.fastq')

#for principle in sampleList:
command = './multiSampleVAFRep.py -i /media/alex/Extra/Dropbox/Code/Sequencing/5.2.2017_InitialAnalysis_4.1.17/2017-05-01_5_0.75 -o /media/alex/Extra/Dropbox/Code/FERMI/varAlleleFrequency -s 4r1.fastq 5r1.fastq 6r1.fastq 7r1.fastq 8r1.fastq 9r1.fastq 10r1.fastq 11r1.fastq 12r1.fastq 13r1.fastq 14r1.fastq 15r1.fastq 16r1.fastq 17r1.fastq 18r1.fastq 20r1.fastq 21r1.fastq 22r1.fastq -p %s -r 0.0015 -c' % ('19r1.fastq')

system(command)
