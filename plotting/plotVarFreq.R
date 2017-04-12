# cool little script that checks if ggplot2 is already installed and if not, installs it
list.of.packages <- c('ggplot2', 'reshape')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dep=TRUE)
library(ggplot2)
library(reshape)

setwd("/media/alex/Extra/Dropbox/Code/R/sequencingWork/4.6.2015 Attempt/")
setwd('/home/alex/Dropbox/Code/R/sequencingWork/4.6.2015 Attempt/')

#subset of the data
#allelefreqshort <- read.table("sortedAlleleFreqs.txt", header = FALSE)
allelefreqshort <- read.table("allelefreqs.txt", header = FALSE)
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

##################
# 3.21.2016 Data #
##################

setwd('/home/alex/Dropbox/Degregori_Lab/5.10.2016/e1r1')
allelefreqshort <- read.table('AF0_plottable.txt', header = TRUE)
hist(allelefreqshort$AO, breaks=100, main='Sample 35 AF=0', xlab='AO', ylab='# Unique Variants')
hist(allelefreqshort$AO, breaks=100, ylim = c(0,10), main='Sample 35 AF=0', xlab='AO', ylab='# Unique Variants')

setwd('/home/alex/Dropbox/Degregori_Lab/5.10.2016/e1r1')
allelefreqshort <- read.table('AF1_plottable.txt', header = TRUE)
hist(allelefreqshort$AO, breaks=100, ylim = c(0,10), main='Sample 35 AF=1', xlab='AO', ylab='# Unique Variants')

setwd('/home/alex/Dropbox/Degregori_Lab/5.10.2016/f1r1')
allelefreqshort <- read.table('AF0_plottable.txt', header = TRUE)
hist(allelefreqshort$AO, breaks=100, main='Sample 36 AF=0', xlab='AO', ylab='# Unique Variants')
hist(allelefreqshort$AO, breaks=100, ylim = c(0,10), main='Sample 36 AF=0', xlab='AO', ylab='# Unique Variants')

setwd('/home/alex/Dropbox/Degregori_Lab/5.10.2016/f1r1')
allelefreqshort <- read.table('AF1_plottable.txt', header = TRUE)
hist(allelefreqshort$AO, breaks=100, ylim = c(0,10), main='Sample 36 AF=1', xlab='AO', ylab='# Unique Variants')

############
# Onc Prob #
############
setwd('~/Desktop')
chanceOnc <- read.table('chanceOnc.txt', header = TRUE)
dat.m <- melt(chanceOnc, 'Sample')
qplot(variable, value, data=dat.m, geom='boxplot',
      main='Freq of Mutation by Region',
      xlab='Genomic Region',
      ylab='Probability of a Base to be Mutated')

####################
# Base Change Bias #
####################
setwd('~/Desktop')
oncoGenes <- read.table('oncoGenes.txt', header = TRUE)

seq=c("A", "C", "T", "G")

df=data.frame('Ref'=rep(seq, each=4), 'Var'=rep(seq, 4), 'Obs'=oncoGenes$Obs)

ggplot(data=df) + aes(x=Var, y=Obs) + 
  geom_bar(stat='identity', position="dodge", fill="grey", color="black")  +
  facet_grid(~Ref) +
  labs(title = "293FT (Samples 33/34) Oncogene Sites")
