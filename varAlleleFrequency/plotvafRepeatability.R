list.of.packages <- c("ggplot2", "ggrepel") # cool little script that checks if ggplot2 is already installed and if not, installs it
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dep=TRUE)

library(ggplot2)
library(ggrepel) # this avoids overlapping labels

# Note that this script cannot begin with a comment for some reason

setwd('/media/alex/Extra/Dropbox/Code/FERMI/varAlleleFrequency')
#setwd('/home/alex/Dropbox/Code/FERMI/varAlleleFrequency')

# This script plots the VAFs of each mutation found between two samples along with a regression
# line and 95% confidence interval in order to understand how repeatable the AFs are for the 
# same variants between samples

oncosites <- c(5073770,7577539,7577119,115256529,115258747,115258744,534287,534288,534289,25398284,25380275,106197266,106197267,106197268,106197269,106155172,106155173,106155174,25457242,25457243,209113112,209113113,90631934,90631838,48649700)

# read in the data
#vafs <- read.table("vafRepeatability.txt", header = TRUE)
vafs <- read.table("outputFile", header = TRUE)
sample1 <- vafs$Sample1
sample2 <- vafs$Sample2
identity <- vafs$Identity

# plot with 95% confidence interval
lm_fit  = lm(sample1 ~ sample2)
x = data.frame(vafs, predict(lm_fit, interval = 'prediction'))

# plot with y=x line
p <- ggplot(vafs, aes(x=sample1, y=sample2, alpha=0.5, label=identity)) +
  geom_point() +
  xlim(0,0.003) +
  ylim(0,0.003) +
  geom_abline(intercept = 0, slope = 1) +
  xlab('Mutation VAFs Sample 503 Deviator Cord Reseq') + ylab('Mutation VAFs Sample 305 Reg Cord Reseq') +
  labs(title = 'Variant Allele Frequencies of Putative Mutations')
# print(p)
jpeg('output2.jpg')
print(p)
dev.off()

# The following produces plot with the oncogegenic mutations identified
'''
library(ggrepel)
vafs$col <- grepl(paste0(oncosites,collapse = "|"), vafs$Identity)
p <- ggplot(vafs, aes(x=Sample1, y=Sample2, alpha=1, color = col)) +
    geom_point() +
    geom_text_repel(aes(label=ifelse(Sample2>0.002 |Sample1>0.002 ,as.character(Identity),''))) +
    xlab('Mutation VAFs 19r1') + ylab('Mutation VAFs Average') +
    labs(title = 'Variant Allele Frequencies of Putative Mutations') 
jpeg('output3.jpg')
print(p)
dev.off()
'''
