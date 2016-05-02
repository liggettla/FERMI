#FERMI
##Fast Extremely Rare Mutation Identification

FERMI is used to identify mutations at an extremely rare frequency.
The program is fully function at the moment, but missing a great deal of 
functionality that will be added in the future to increase useability.

####Usage Instructions
Please see the [FERMI Homepage](http://liggettla.github.io/FERMI/) for
download and usage instructrions. This page will be progressively updated
to be as complete as possible.






*"Where is everybody?!"*

####Help File
optional arguments:
  -h, --help            show this help message and exit
  --nfo NFO, -n NFO     Info writeup about a particular run that will be
                        output in the run directory.
  --largefiles, -l      Outputs all generated fastq files generated during
                        analysis.
  --avoidalign, -a      Only runs through initial analysis of input fastq
                        files, and does not align to reference or call
                        variants.
  --outdir OUTDIR, -o OUTDIR
                        Specifies output directory where all analysis files
                        will be written
  --indir INDIR, -i INDIR
                        Specifies the input directory that contains the fastq
                        files to be analyzed.
  --single, -s          Only process a single set of paired end reads.
  --prevdict PREVDICT, -p PREVDICT
                        Specify a previously output pickle file containing
                        collapsed fastq data as an input instead of raw fastq
                        files.
  --umimismatch UMIMISMATCH, -u UMIMISMATCH
                        Specify the number of mismatches allowed in a UMI pair
                        to still consider as the same UMI
  --varthresh VARTHRESH, -v VARTHRESH
                        Specify the percentage of reads that must contain a
                        particular base for that base to be used in the final
                        consensus read
  --readsupport READSUPPORT, -r READSUPPORT
                        Specifies the number of reads that must have a given
                        UMI sequence in order to be binned as a true capture
                        event, and not be thrown out.
  --clustersubmit, -c   Submit run to cluster computing rather than running
                        locally
