#BSUB -J FERMI[1]
#BSUB -R "span[hosts=1]"
#BSUB -n 1
#BSUB -R "rusage[mem=44445]"
#BSUB -o /vol3/home/liggettl/logs/output.%J -e /vol3/home/liggettl/logs/error.%J job01

python main.py
