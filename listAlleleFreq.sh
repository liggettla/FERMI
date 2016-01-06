#the purpose of this script is to genereate a file containing allele freqs
#to be used in plotting

files=file

VCFDIR=/dir

currentfile=$files

#looks for lines starting 'chr' and then grabs the values in columns
#1,2,4, and 5 then searches for RO and AO and calculates frequency of
#allele variant by AO/(AO+RO)
awk 'v {
        ro=gensub(/^.*;RO=([0-9]*).*$/, "\\1", "1");
        printf "%s %f\n", f, (ao/(ao + ro)); v=0
    }
    /^chr/ {ao=gensub(/^.*;AO=([0-9]*).*$/,"\\1", "1");
        v=1;
        f=$1 FS $2 FS $4 FS $5
        }' $VCFDIR/$currentfile.vcf >> $VCFDIR/totalvariants.txt #output all to single file
#        }' $VCFDIR/$currentfile.vcf > $VCFDIR/$currentfile.sorted.txt #output to indiv files

awk '{print $5}' $VCFDIR/totalvariants.txt > $VCFDIR/allelefreqs.txt
