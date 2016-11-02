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

# read in the data
vafs <- read.table("vafRepeatability.txt", header = TRUE)
#vafs <- read.table("outputFile", header = TRUE)
sample1 <- vafs$Sample1
sample2 <- vafs$Sample2
identity <- vafs$Identity

# plot with 95% confidence interval
lm_fit  = lm(sample1 ~ sample2)
x = data.frame(vafs, predict(lm_fit, interval = 'prediction'))

# this appears to do the same as the above
# coef(lm(sample1 ~ sample2, data = vafs))
# geom_smooth(method = "lm", se = TRUE) +

p <- ggplot(x, aes(x=sample1, y=sample2, alpha=0.5, label=identity)) +
  #geom_text(aes(label=identity),hjust='inward', vjust='inward', angle=0) + # this labels all points
  #geom_text(aes(label=ifelse(sample2>0.002|sample1>0.002,as.character(identity),'')),hjust=0,vjust=0) + # this labels points above particular frequency
  #geom_text_repel(aes(label=ifelse(sample2>0.002|sample1>0.002,as.character(identity),''))) + # this labels points above freq and does not allow overlap
  xlab('Mutation VAFs C1 F34') +
  ylab('Muation VAFs E1 F41') +
  labs(title = 'Variant Allele Frequencies of Putative Mutations (W/ Regression)') +
  geom_point() +
  geom_smooth(method = 'lm', aes(fill = 'confidence'), alpha = 0.15) +
  scale_fill_manual('Interval', values = c('green', 'blue'))
# print(p)
jpeg('output1.jpg')
print(p)
dev.off()

p <- ggplot(vafs, aes(x=sample1, y=sample2, alpha=0.5, label=identity)) +
  geom_point() +
  #geom_text(aes(label=identity), hjust='inward', vjust='inward', angle=0) + # this labels all points
  #geom_text(aes(label=ifelse(sample2>0.002|sample1>0.002,as.character(identity),'')),hjust=0,vjust=0) + # this labels points above particular frequency
  geom_text_repel(aes(label=ifelse(sample2>0.4|sample1>0.4,as.character(identity),''))) + # this labels points above freq and does not allow overlap
  geom_abline(intercept = 0, slope = 1) +
  xlab('Mutation VAFs C1 F34') +
  ylab('Muation VAFs E1 F41') +
  labs(title = 'Variant Allele Frequencies of Putative Mutations')
# print(p)
jpeg('output2.jpg')
print(p)
dev.off()
