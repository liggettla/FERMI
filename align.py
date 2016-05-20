####################
#Align to Reference#
####################
# MAKE SURE TO FIRST RUN THE FOLLOWING COMMAND:
# module load bwa

# @RG for readgroup ID
# -R proceed with suboptimal aligments, improves paired end mapping at the
# cost of speed
# -U 40 increases penalty for unpaired readpair b/c ignoring put right after -M
# translocations, just looking for single variants
# -M outputs for picard compatibility
# -t 13 uses 13 threads (24 per node so uses separate node per sample)
# -bS - this reads in the binary format sam file and converts to bam file
# dash tells standard in to read directly from the stream

# bwa mem -R '@RG\tID:'$currentfile.fastq'\tSM:'$currentfile.fastq \
#    -M $REF $R1 $R2 -t 13 \
#    | samtools view -bS - \
#    | samtools sort - $RESDIR/$currentfile.fastq

# only run this on the cluster

from os import system

def align(final_output_file, REF):
    bamOut = final_output_file.strip('fastq') + 'bam'

    # this all-in-one call seems to be problematic with the new cluster memory limitations
    #system("bwa mem %s %s | samtools view -bS - | samtools sort > %s" % (REF, final_output_file, bamOut))
    system("bwa mem %s %s > temp.sam" % (REF, final_output_file))
    system("samtools view -bS temp.sam > temp.txt")
    system("samtools sort temp.txt > %s" % (bamOut))
    system("samtools index %s" % (bamOut))
    system("rm temp.sam")
    system("rm temp.txt")

    return bamOut
