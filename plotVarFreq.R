#library(ggplot2)

setwd("/media/alex/Extra/Dropbox/Code/R/Sequencing Work/4.6.2015 Attempt/")

#subset of the data
allelefreqshort <- read.table("sortedAlleleFreqs.txt", header = FALSE)
hist(log10(allelefreqshort$V1), xlim = c(-4,0), breaks=100)

hist(log(allelefreqshort$AlleleFreq), xlim = c(-15,0), breaks=100)

#all of the data
allelefreq <- read.table("allelefreqs.txt", header = TRUE)
hist(log10(allelefreq$AlleleFreq), xlim = c(-2,0), breaks=100)
hist(log10(allelefreq$AlleleFreq), xlim = c(-1,0), breaks=300)

#HiSeq data
#Each different set is using different cut-offs for the number of times a variant is present in a given capture
setwd("/media/alex/Extra/Dropbox/TrueSeqPipeline/results/5_Reads_75_Percent")
allelefreqshort <- read.table("allelefreqs.txt", header = FALSE)
finalPlot = hist(log10(allelefreqshort$V1), xlim = c(-4,0), breaks=100, main='5 Reads 75 Percent', xlab='Variant Frequency', ylab='Number of Variants')

setwd("/media/alex/Extra/Dropbox/TrueSeqPipeline/results/10_Reads_95_Percent")
allelefreqshort <- read.table("allelefreqs.txt", header = FALSE)
hist(log10(allelefreqshort$V1), xlim = c(-4,0), breaks=100, main='10 Reads 95 Percent', xlab='Variant Frequency', ylab='Number of Variants')

setwd("/media/alex/Extra/Dropbox/TrueSeqPipeline/results/10_Reads_99_Percent")
allelefreqshort <- read.table("allelefreqs.txt", header = FALSE)
hist(log10(allelefreqshort$V1), xlim = c(-0.4,0), breaks=100, main='10 Reads 99 Percent', xlab='Variant Frequency', ylab='Number of Variants')
