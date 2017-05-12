list.of.packages <- c("ggplot2", "ggrepel") # cool little script that checks if ggplot2 is already installed and if not, installs it
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dep=TRUE)

library(ggplot2)
library(ggrepel) # this avoids overlapping labels

# Note that this script cannot begin with a comment for some reason

setwd('/media/alex/Extra/Dropbox/Code/FERMI/varAlleleFrequency')
#setwd('/media/alex/Extra/Dropbox/Code/FERMI/varAlleleFrequency')
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

# this appears to do the same as the above
# coef(lm(sample1 ~ sample2, data = vafs))
# geom_smooth(method = "lm", se = TRUE) +

# plot with regression line
#p <- ggplot(x, aes(x=sample1, y=sample2, alpha=0.5, label=identity)) +
    #geom_text(aes(label=identity),hjust='inward', vjust='inward', angle=0) + # this labels all points
    #geom_text(aes(label=ifelse(sample2>0.002|sample1>0.002,as.character(identity),'')),hjust=0,vjust=0) + # this labels points above particular frequency
    #geom_text_repel(aes(label=ifelse(sample2>0.002|sample1>0.002,as.character(identity),''))) + # this labels points above freq and does not allow overlap
#  xlab('Mutation VAFs C1 F34') +
#  ylab('Muation VAFs E1 F41') +
#  labs(title = 'Variant Allele Frequencies of Putative Mutations (W/ Regression)') +
#  geom_point() +
#  geom_smooth(method = 'lm', aes(fill = 'confidence'), alpha = 0.15) +
#  scale_fill_manual('Interval', values = c('green', 'blue'))
# print(p)
#jpeg('output1.jpg')
#print(p)
#dev.off()

# the following line finds the maximum vafs for use as plotting cutoffs
# awk -v max1=0.0 -v max2=0.0 '{if($1!="Sample1" && $1>max1 && 2>max2){max1=$1; max2=$2}}END{print max1; print max2}' vafRepeatability.txt
# plot with y=x line
p <- ggplot(vafs, aes(x=sample1, y=sample2, alpha=0.5, label=identity)) +
  geom_point() +
  #xlim(0,1) +
  #ylim(0,1) +
  #geom_text(aes(label=identity), hjust='inward', vjust='inward', angle=0) + # this labels all points
  #geom_text(aes(label=ifelse(sample2>0.002|sample1>0.002,as.character(identity),'')),hjust=0,vjust=0) + # this labels points above particular frequency
  
    #geom_text_repel(aes(label=ifelse(sample2>0.4|sample1>0.4,as.character(identity),''))) + # this labels points above freq and does not allow overlap
    #geom_text_repel(aes(label=ifelse(sample2>0.002|sample1>0.002,as.character(identity),''))) + # this labels points above freq and does not allow overlap
    geom_text_repel(aes(label=ifelse(sample2>0.007 |sample1>0.007 ,as.character(identity),''))) + # this labels points above freq and does not allow overlap
  geom_abline(intercept = 0, slope = 1) +
    xlab('Mutation VAFs 23r1') + ylab('Mutation VAFs 2r1') +
    #xlab('Mutation VAFs A1 305 Cord') + ylab('Muation VAFs B1 305 Cord') +
    #xlab('Mutation VAFs C1 300 F34') + ylab('Muation VAFs D1 300 F34') +
    #xlab('Mutation VAFs E1 301 F41') + ylab('Muation VAFs F1 301 F41') +
    #xlab('Mutation VAFs G1 302 F34') + ylab('Muation VAFs H1 302 F34') +
    #xlab('Mutation VAFs A2 303 F46') + ylab('Muation VAFs B2 303 F46') +
    #xlab('Mutation VAFs C2 304 M30') + ylab('Muation VAFs D2 304 M30') +
  labs(title = 'Variant Allele Frequencies of Putative Mutations')
# print(p)
jpeg('output2.jpg')
print(p)
dev.off()

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
