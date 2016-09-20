#BSUB -J MERGE 
#BSUB -R "span[hosts=1]"
#BSUB -n 1
#BSUB -R "rusage[mem=10]"
#BSUB -o /vol3/home/liggettl/logs/output.%J -e /vol3/home/liggettl/logs/error.%J job01

#This is a script to concatenate and extract fastq files in batch

cd "/vol3/home/liggettl/sequencingData/8.3.2016" || exit 1

# example filename: A1_S1_L001_R1_001.fastq.gz
for num in {1..3}; do
    for letter in {A..H}; do
        for subgroup in R1 R2; do
            zcat "$letter$num"_*_L*_"$subgroup"_001.fastq.gz > "$letter$num-$subgroup".txt
        done
    done
done
