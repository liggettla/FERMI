#BSUB -J FERMI[1]
#BSUB -e ~/logs/FERMI%I.%J.err
#BSUB -o ~/logs/FERMI%I.%J.out
#BSUB -R "span[hosts=1]"
#BSUB -n 1

python main.py
