#BSUB -J UMICon[1]
#BSUB -e ~/logs/dicCol.%I.%J.err
#BSUB -o ~/logs/dicCol.%I.%J.out
#BSUB -R "span[hosts=1]"
#BSUB -n 1


python main.py
