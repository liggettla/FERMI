#FERMI
##Fast Extremely Rare Mutation Identification

FERMI is used to identify mutations at an extremely rare frequency.
The program is fully function at the moment, but missing a great deal of 
functionality that will be added in the future to increase useability.

####Usage Instructions
Please see the [FERMI Homepage](http://liggettla.github.io/FERMI/) for
download and usage instructrions. This page will be progressively updated
to be as complete as possible.

##Installation
The following instructions are for Debian Linux, but with some adjustment should work anywhere.

Download FERMI
`git clone https://github.com/liggettla/FERMI`

Install bwa
`sudo apt-get install bwa`

Download freebayes
`git clone --recursive git://github.com/ekg/freebayes.git`
`cd freebayes`
`make`

##Example Runs
Take input files from testInput and put analyzed files in testOutput:

`./fermi.py -i testInput/ -o testOutput/`

Submit run to computing cluster:

`./fermi.py -c -i testInput/ -o testOutput/`

Allow 2 UMI mismatches, 0.6 variant threshold and 3 supporting reads of a UMI pair:

`./fermi.py -u 2 -v 0.6 -r 3 -i testInput/ -o testOutput/`

Output all fastq files generated in the analysis:

`./fermi.py -l -i testInput/ -o testOutput/`

Avoid alignment and variant calling for testing:

`./fermi.py -a -i testInput/ -o testOutput/`

*"Where is everybody?!"*

####Help File
optional arguments:

  &nbsp;&nbsp;&nbsp;&nbsp;-h, --help            show this help message and exit
  
  &nbsp;&nbsp;&nbsp;&nbsp;--nfo NFO, -n NFO     Info writeup about a particular run that will be
                        output in the run directory.
                        
 &nbsp;&nbsp;&nbsp;&nbsp; --largefiles, -l      Outputs all generated fastq files generated during
                        analysis.
                        
 &nbsp;&nbsp;&nbsp;&nbsp; --avoidalign, -a      Only runs through initial analysis of input fastq
                        files, and does not align to reference or call
                        variants.
                        
  &nbsp;&nbsp;&nbsp;&nbsp;--outdir OUTDIR, -o OUTDIR
                        Specifies output directory where all analysis files
                        will be written
                        
  &nbsp;&nbsp;&nbsp;&nbsp;--indir INDIR, -i INDIR
                        Specifies the input directory that contains the fastq
                        files to be analyzed.
                        
  &nbsp;&nbsp;&nbsp;&nbsp;--single, -s          Only process a single set of paired end reads.
  
  &nbsp;&nbsp;&nbsp;&nbsp;--prevdict PREVDICT, -p PREVDICT
                        Specify a previously output pickle file containing
                        collapsed fastq data as an input instead of raw fastq
                        files.
                        
  &nbsp;&nbsp;&nbsp;&nbsp;--umimismatch UMIMISMATCH, -u UMIMISMATCH
                        Specify the number of mismatches allowed in a UMI pair
                        to still consider as the same UMI
                        
  &nbsp;&nbsp;&nbsp;&nbsp;--varthresh VARTHRESH, -v VARTHRESH
                        Specify the percentage of reads that must contain a
                        particular base for that base to be used in the final
                        consensus read
                        
  &nbsp;&nbsp;&nbsp;&nbsp;--readsupport READSUPPORT, -r READSUPPORT
                        Specifies the number of reads that must have a given
                        UMI sequence in order to be binned as a true capture
                        event, and not be thrown out.
                        
  &nbsp;&nbsp;&nbsp;&nbsp;--clustersubmit, -c   Submit run to cluster computing rather than running
                        locally
