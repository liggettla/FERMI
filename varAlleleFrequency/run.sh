#!/bin/bash

#./multiSampleVAFRep.py -i /home/alex/Dropbox/Degregori_Lab/Experiments/2.13.2017_293_Comparisons/2017-01-09_5_0.55/ -o /home/alex/Dropbox/Code/FERMI/varAlleleFrequency -s e1r1.fastq f1r1.fastq -p d1c1r1.fastq -r 0.003
#./multiSampleVAFRep.py -i /home/alex/Dropbox/Degregori_Lab/Experiments/2.13.2017_293_Comparisons/2017-01-09_5_0.55 -o /home/alex/Dropbox/Code/FERMI/varAlleleFrequency -s e1r1.fastq -p f1r1.fastq -r 0.003

./multiSampleVAFRep.py -i /home/alex/Dropbox/Degregori_Lab/Experiments/2.13.2017_293_Comparisons/compareSubstitution/ -o /home/alex/Dropbox/Code/FERMI/varAlleleFrequency/ -s subst -p noSubst -r 0.5

#./martincorena.py -i /media/alex/Extra/Dropbox/Degregori_Lab/Experiments/10.18.2016_Martincorena -o /media/alex/Extra/Dropbox/Code/FERMI/varAlleleFrequency -s Sample1 Sample4 Sample2 -p Sample3
