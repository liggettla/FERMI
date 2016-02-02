R1.fastq and R2.fastq have five unique captures in them.
The first two are repeated with a different UMI sequence in R1.fastq while there is only a single unique
capture of the third sequence. There are 10 repetitions of each read.
This should allow each of the captures to pass the default fileters but the third read should hopefully be
automatically ignored eventually because there was only a single unique capture.
