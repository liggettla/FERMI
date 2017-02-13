#BSUB -J counts
#BSUB -R "span[hosts=1]"
#BSUB -n 1
#BSUB -R "rusage[mem=30]"
#BSUB -o /vol3/home/liggettl/logs/output.%J -e /vol3/home/liggettl/logs/error.%J job01

python countStrands.py -i /vol3/home/liggettl/sequencingData/8.3.2016/unzipped/
