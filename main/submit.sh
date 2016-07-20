#BSUB -J FERMI[1]
#BSUB -e ~/logs/FERMI%I.%J.err
#BSUB -o ~/logs/FERMI%I.%J.out
#BSUB -R "span[hosts=1]"
#BSUB -n 1
#BSUB -R "rusage[mem=55]"

#########
# fermi #
#########
#module load vt
#module load bwa

./fermi.py -i ./testInput -o ./testOutput -u 1 -v 0.51 -r 3 -f 1 -n 'This is some information.'
#./fermi.py -i ./testInput -o ./testOutput -u 1 -v 0.51 -r 3 -c -f 1
#./fermi.py -i /vol3/home/liggettl/sequencingData/3.21.2016/unzipped/untrimmed/unzipped -o /vol3/home/liggettl/sequencingData/3.21.2016/analysis -u 1 -v 0.51 -r 3 -c

#~/FERMI/fermi.py -i /vol3/home/liggettl/sequencingData/3.21.2016/unzipped/untrimmed/unzipped -o /vol3/home/liggettl/sequencingData/3.21.2016/analysis -u 1 -v 0.51 -r 3 -c -f 1

########################
# sample repeatability #
########################
#./sampleRepeatability.py -i /vol3/home/liggettl/sequencingData/3.21.2016/analysis/2016-05- 03_3_0.51 -o /vol3/home/liggettl/sequencingData/3.21.2016/analysis/2016-05-03_3_0.51/BlockDecomposed -s a1r1.fastq b1r1.fastq

#############
# base bias #
#############
#./baseChangeBias.py -i /vol3/home/liggettl/sequencingData/3.21.2016/analysis/2016-05-03_3_0.51 -o /vol3/home/liggettl/sequencingData/3.21.2016/analysis/2016-05-03_3_0.51 -s f1r1.fastq