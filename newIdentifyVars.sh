currentfile=finalOutput.bam
VCFDIR=/vol3/home/liggettl/TruSeqPanel/8.17.2015_HiSeqFastq/combinedUMI/vcfs

awk '/^chr/{match($0,/;AO=([^;]*);.*;RO=([^;]*);/,v); 
            print $1,$2,$4,$5,v[1]/(v[1]+v[2])}' $VCFDIR/$currentfile.vcf > $VCFDIR/totalvariants.txt #output all to single file

awk '{print $5}' $VCFDIR/totalvariants.txt > $VCFDIR/allelefreqs.txt
