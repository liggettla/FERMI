list.of.packages <- c("ggplot2", "ggrepel") # cool little script that checks if ggplot2 is already installed and if not, installs it
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dep=TRUE)

library(ggplot2)
library(ggrepel) # this avoids overlapping labels

# Note that this script cannot begin with a comment for some reason

#setwd('/media/alex/Extra/Dropbox/Code/FERMI/varAlleleFrequency')
setwd('/home/alex/Dropbox/Code/FERMI/varAlleleFrequency')

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
p <- ggplot(vafs, aes(x=sample1, y=sample2, alpha=0.5, label=identity, size=15)) +
  geom_point() +
  xlim(0,0.003) +
  ylim(0,0.003) +
  #xlim(0,1) +
  #ylim(0,1) +
  geom_abline(intercept = 0, slope = 1, size=3)+ # y=x line
  xlab('\nVAF Individual 15') + ylab('VAF Mean\n') +
  #labs(title = 'C-T Variants Coding Regions\n')+
# 2a
  #xlab('\nVAF Individual 15') + ylab('VAF Individual 7\n') +
# 2b
  xlab('\nVAF Individual 15') + ylab('VAF Mean\n') +
# 2.2.3
  #xlab('\nVAF Individual 7') + ylab('VAF Individual 7\n') +
# 3a
  #xlab('\nVAF HCT116 MMR-') + ylab('VAF HCT116 MMR+\n') +
# 3b
  #xlab('\nVAF HCT116 MMR-') + ylab('VAF Mean\n') +
# 3c
  #xlab('\nVAF Individual 19') + ylab('VAF Mean\n') +
# 3d
  #xlab('\nVAF Individual 2') + ylab('VAF Mean\n') +
# 3e
  #xlab('\nVAF Individual 2') + ylab('VAF Individual 19\n') +
# 3f
  #xlab('\nVAF HCT116 MMR-') + ylab('VAF Mean Individuals 19, 2\n') +
# 3h
  #xlab('\nVAF HCT116 MMR-') + ylab('VAF Individual 2\n') +
  #labs(title = 'C-N/G-N Variants')+
# 3i
  #xlab('\nVAF HCT116 MMR-') + ylab('VAF Individual 2\n') +
  #labs(title = 'T-N/A-N Variants')+
# 2s1
  #xlab('\nVAF Individual 15') + ylab('VAF Individual 7\n') +
# 3s2
  #xlab('\nVAF Individual 15') + ylab('VAF Individual 7\n') +
  #labs(title = 'Total Variants\n')+
  #labs(title = 'C-T/G-A Variants\n')+
  #labs(title = 'T-C/A-G Variants\n')+
  #labs(title = 'C-A/G-T Variants\n')+
  #labs(title = 'T-A/A-T Variants\n')+
  #labs(title = 'C-G/G-C Variants\n')+
  #labs(title = 'T-G/A-C Variants\n')+
# 3s3
  #xlab('\nVAF HCT116 MMR-') + ylab('VAF HCT116 MMR+\n') +
  #labs(title = 'Total Variants\n')+
  #labs(title = 'C-T/G-A Variants\n')+
  #labs(title = 'C-G/G-C Variants\n')+
  #labs(title = 'C-A/G-T Variants\n')+
  #labs(title = 'T-C/A-G Variants\n')+
  #labs(title = 'T-G/A-C Variants\n')+
  #labs(title = 'T-A/A-T Variants\n')+
# 3s4
  #xlab('\nVAF Individual 2') + ylab('VAF Mean\n') +
  #labs(title = 'Total Variants\n')+
  #labs(title = 'C-T/G-A Variants\n')+
  #labs(title = 'C-G/G-C Variants\n')+
  #labs(title = 'C-A/G-T Variants\n')+
  #labs(title = 'T-C/A-G Variants\n')+
  #labs(title = 'T-G/A-C Variants\n')+
  #labs(title = 'T-A/A-T Variants\n')+
# 3s5
  #xlab('\nVAF Individual 19') + ylab('VAF Mean\n') +
  #labs(title = 'Total Variants\n')+
  #labs(title = 'C-T/G-A Variants\n')+
  #labs(title = 'C-G/G-C Variants\n')+
  #labs(title = 'C-A/G-T Variants\n')+
  #labs(title = 'T-C/A-G Variants\n')+
  #labs(title = 'T-G/A-C Variants\n')+
  #labs(title = 'T-A/A-T Variants\n')+

  #geom_smooth(method=lm, se=FALSE, size=3)+ # regression line
  theme_bw()+ # no gray background
  theme(panel.border = element_blank())+ # no border
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ # no gridlines
  theme(axis.title = element_text(size = 50))+ # change label size
  theme(plot.title = element_text(size = 50))+ # change title size
  theme(plot.title = element_text(hjust = 0.5))+ # center title
  theme(axis.text.x = element_text(size = 50, colour="black", angle=90))+ # change tick size
  theme(axis.text.y = element_text(size = 50, colour="black"))+ # change tick size
  theme(legend.position="none")+ # no legend
  theme(axis.ticks = element_line(colour = "black", size = 2))+ # hide ticks
  theme(axis.line = element_line(colour = "black", size=3)) # add axis
# print(p)
#jpeg("output2.jpg")
jpeg("output2.jpg", units="in", width=17, height=17, res=500)
print(p)
dev.off()

'''
m <- lm(sample2 ~ sample1);
print(m)

# The following produces plot with the oncogegenic mutations identified
library(ggrepel)
vafs$col <- grepl(paste0(oncosites,collapse = "|"), vafs$Identity)
p <- ggplot(vafs, aes(x=Sample1, y=Sample2, alpha=1, color = col)) +
    geom_point(size=5) +
    geom_abline(intercept = 0, slope = 1, size=3)+ # y=x line
    xlab("\nVAF Individual 2") + ylab("VAF Mean\n") +
    xlim(0,0.006) +
    ylim(0,0.006) +
    #geom_text_repel(aes(label=ifelse(Sample2>0.002 |Sample1>0.002 ,as.character(Identity),""))) +
    theme_bw()+ # no gray background
    theme(panel.border = element_blank())+ # no border
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ # no gridlines
    theme(axis.title = element_text(size = 50))+ # change label size
    theme(plot.title = element_text(size = 50))+ # change title size
    theme(plot.title = element_text(hjust = 0.5))+ # center title
    theme(axis.text.x = element_text(size = 50, colour="black", angle=90))+ # change tick size
    theme(axis.text.y = element_text(size = 50, colour="black"))+ # change tick size
    theme(legend.position="none")+ # no legend
    theme(axis.ticks = element_line(colour = "black", size = 2))+ # hide ticks
    theme(axis.line = element_line(colour = "black", size=3)) # add axis
jpeg("output2.jpg", units="in", width=17, height=17, res=500)
print(p)
dev.off()
'''
