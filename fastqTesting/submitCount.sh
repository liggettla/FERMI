#BSUB -J COUNT
#BSUB -R "span[hosts=1]"
#BSUB -n 1
#BSUB -R "rusage[mem=5]"
#BSUB -o /vol3/home/liggettl/logs/output.%J -e /vol3/home/liggettl/logs/error.%J job01

#submit with: bsub -n 1 < submitCount.sh
python countACTG.py -i /vol3/home/liggettl/sequencingData/3.21.2016/unzipped/untrimmed/unzipped/ -o /vol3/home/liggettl/FERMI/fastqTesting/ -f e1r1.fastq 
