#BSUB -J FERMI[1]
#BSUB -e ~/logs/FERMI%I.%J.err
#BSUB -o ~/logs/FERMI%I.%J.out
#BSUB -R "span[hosts=1]"
#BSUB -n 1
#BSUB -R "rusage[mem=40]"

#########
# fermi #
#########
#module load vt
#module load bwa

###############
# Run Locally #
###############
# this will avoid aligning
# ./fermi.py -i /media/alex/Extra/Dropbox/Code/FERMI/testInput -o /media/alex/Extra/Dropbox/Code/FERMI/testOutput -q 50 -e -u 1 -v 0.8 -r 5 -f 5 -a -n 'This is some information.'

# this will run a complete run
./fermi.py -i /media/alex/Extra/Dropbox/Code/FERMI/testInput -o /media/alex/Extra/Dropbox/Code/FERMI/testOutput -q 120 -ex -u 1 -v 0.55 -r 5 -f 5 -b 'freebayes' -y '/media/alex/Extra/Dropbox/Code/ReferenceGenomes/hg19.fa' -n 'This is some information.'

###############
# Cluster Run #
###############
#~/FERMI/main/fermi.py -i /vol3/home/liggettl/sequencingData/8.3.2016/unzipped -o /vol3/home/liggettl/sequencingData/8.3.2016/analysis -xce -u 1 -v 0.95 -r 5 -f 5 -q 120 -y '/vol3/home/liggettl/refgenomes/hg19.fa' -n 'The purpose of this run is to compare variant dropout rates with previous data but instead of using the artificial hg19, this is using the original hg19.'

#################
# 8.3.2016 data #
#################
#~/FERMI/main/fermi.py -i /vol3/home/liggettl/sequencingData/8.3.2016/unzipped -o /vol3/home/liggettl/sequencingData/8.3.2016/analysis -u 1 -v 0.6 -r 3 -c -f 1 -e -n 'This is the first run of the 8.3.2016 dataset. Everything is being run with routine settings for an initial look at the data.'

########################
# sample repeatability #
########################
#./sampleRepeatability.py -i /vol3/home/liggettl/sequencingData/3.21.2016/analysis/2016-05- 03_3_0.51 -o /vol3/home/liggettl/sequencingData/3.21.2016/analysis/2016-05-03_3_0.51/BlockDecomposed -s a1r1.fastq b1r1.fastq

#############
# base bias #
#############
#./baseChangeBias.py -i /vol3/home/liggettl/sequencingData/3.21.2016/analysis/2016-05-03_3_0.51 -o /vol3/home/liggettl/sequencingData/3.21.2016/analysis/2016-05-03_3_0.51 -s f1r1.fastq
