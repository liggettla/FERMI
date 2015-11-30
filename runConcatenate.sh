#BSUB -J UMICon[1]
#BSUB -e ~/logs/dicCol.%I.%J.err
#BSUB -o ~/logs/dicCol.%I.%J.out
#BSUB -R "span[hosts=1]"
#BSUB -n 1

##############
#Lab Computer#
##############
python /media/alex/Extra/Dropbox/Code/TrueSeqMiSeqScripts/hiSeqAnalysis/main.py

#####################
#Cluster Single File#
#####################
#python /vol3/home/liggettl/TruSeqPanel/Scripts/hiSeqAnalysis/main.py
