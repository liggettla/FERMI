#BSUB -J align[1]
#BSUB -e logs/align.%I.%J.err
#BSUB -o logs/align.%I.%J.out
#BSUB -R "span[hosts=1]"
#BSUB -n 1

#Reference genome must be indexed for bwa to run
#Do this with bwa index hg19.fa then refer to hg19.fa in the bwa mem line
#that has all of the index files with it in the same directory 
#bwa help here: http://bio-bwa.sourceforge.net/bwa.shtml

set -exo pipefail

#load these on headnode if running on LRS cluster
#module load bwa
#module load java/1.7

#this is the location of the fastqs, picard, and ref genome
#FQDIR=/vol3/home/liggettl/TruSeqPanel/8.17.2015_HiSeqFastq/25_Reads_75_Percent
#RESDIR=/vol3/home/liggettl/TruSeqPanel/8.17.2015_HiSeqFastq/25_Reads_75_Percent
PICARD=/vol3/home/liggettl/ExomeSeq/bin/picard-tools-1.83
#REF=/vol3/home/liggettl/refgenomes/hg19.fa
REF=/media/alex/Extra/Dropbox/Code/ReferenceGenomes/hg19.fa
#BAMDIR=/vol3/home/liggettl/TruSeqPanel/8.17.2015_HiSeqFastq/25_Reads_75_Percent
#VCFDIR=/vol3/home/liggettl/TruSeqPanel/8.17.2015_HiSeqFastq/25_Reads_75_Percent

#this assigns the passed output folder variable to FQDIR
#in the format of /outputDir
FQDIR=$1
RESDIR=$FQDIR
BAMDIR=$FQDIR
VCFDIR=$FQDIR

#mkdir -p $FQDIR

#currentfile=${files[$(($LSB_JOBINDEX - 1))]}
currentfile=finalOutput
R1=$FQDIR/$currentfile.fastq

#@RG for readgroup ID
#-R proceed with suboptimal aligments, improves paired end mapping at the
#cost of speed
#-U 40 increases penalty for unpaired readpair b/c ignoring put right after -M
#translocations, just looking for single variants
#-M outputs for picard compatibility
#-t 13 uses 13 threads (24 per node so uses separate node per sample)
#-bS - this reads in the binary format sam file and converts to bam file
#dash tells standard in to read directly from the stream

#bwa mem -R '@RG\tID:'$currentfile.fastq'\tSM:'$currentfile.fastq \
#    -M $REF $R1 $R2 -t 13 \
#    | samtools view -bS - \
#    | samtools sort - $RESDIR/$currentfile.fastq

#For Paired-end Reads
#bwa mem $REF $R1 $R2 | samtools view -bS - | samtools sort - $RESDIR/$currentfile.bam

###########
#Alignment#
###########
#bam filenames
files=finalOutput
#files=(UMIoutput_sample5_goodCollapse.bam
#UMIoutput_sample8_R2_goodCollapse.bam)

#currentfile=${files[$(($LSB_JOBINDEX - 1))]}
currentfile=$files

#For Single-end Reads
#bwa mem is better than bwa for reads > 70bp
bwa mem $REF $R1 | samtools view -bS - | samtools sort - $RESDIR/$currentfile

samtools index $RESDIR/$currentfile.bam

###############
#Call Variants#
###############

#currentfile=${files[$(($LSB_JOBINDEX - 1))]}
currentfile=$files

#calls variants only if they are at least 0.00001% of the calls 1/10^6
#this should thus output every variant found in the alignment
#by default freebayes uses 0.2 (20%)

#Cluster
#/vol3/home/liggettl/TruSeqPanel/Scripts/freebayes/freebayes -F 0.0000001 --fasta-reference $REF $BAMDIR/$currentfile.bam > $VCFDIR/$currentfile.vcf

#Lab
/media/alex/Extra/Dropbox/Code/freebayes -F 0.0000001 --fasta-reference $REF $BAMDIR/$currentfile.bam > $VCFDIR/$currentfile.vcf
