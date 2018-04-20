list.of.packages <- c("ggplot2", "ggrepel") # cool little script that checks if ggplot2 is already installed and if not, installs it
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dep=TRUE)

library(ggplot2)
library(ggrepel) # this avoids overlapping labels

# Note that this script cannot begin with a comment for some reason

setwd('/media/alex/Extra/Dropbox/Code/FERMI/varAlleleFrequency')
#setwd('/media/alex/mainstorage/LinuxDropbox/Dropbox/Code/FERMI/varAlleleFrequency')

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

lm_eqn <- function(df,y,x){
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(coef(m)[1], digits = 2), 
              b = format(coef(m)[2], digits = 2), 
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}

# plot with y=x line
p <- ggplot(vafs, aes(x=sample1, y=sample2, alpha=0.5, label=identity, size=15)) +
  geom_point() +
  geom_text(aes(x = 0.0015, y = 0.003, label = lm_eqn(vafs, Sample1, Sample2)), parse = TRUE)+ # regression formula
  xlim(0,0.003) +
  ylim(0,0.003) +
  xlim(0,0.01) +
  ylim(0,0.01) +
  xlim(0,1) +
  ylim(0,1) +
  geom_abline(intercept = 0, slope = 1, size=3)+ # y=x line
  #xlab('\nVAF Individual 15') + ylab('VAF Mean\n') +
  #labs(title = 'C-T Variants Coding Regions\n')+
# 2a
  #xlab('\nVAF Individual 15') + ylab('VAF Individual 7\n') +
# 2b
  #xlab('\nVAF Individual 15') + ylab('VAF Mean\n') +
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
  #xlab('\nVAF Individual 2') + ylab('VAF Mean\n') +
  #labs(title = 'C>N/G>N Variants')+
# 3i
  #xlab('\nVAF Individual 2') + ylab('VAF Mean\n') +
  #labs(title = 'T>N/A>N Variants')+
# 2s1
  #xlab('\nVAF Individual 15') + ylab('VAF Individual 7\n') +
# s10 
  #xlab('\nVAF Individual 15') + ylab('VAF Individual 7\n') +
  #labs(title = 'Total Variants\n')+
  #labs(title = 'C-T/G-A Variants\n')+
  #labs(title = 'T-C/A-G Variants\n')+
  #labs(title = 'C-A/G-T Variants\n')+
  #labs(title = 'T-A/A-T Variants\n')+
  #labs(title = 'C-G/G-C Variants\n')+
  #labs(title = 'T-G/A-C Variants\n')+
# s4
  #xlab('\nAllele 1') + ylab('Allele 2\n') +
# s11 
  #xlab('\nVAF HCT116 MMR-') + ylab('VAF HCT116 MMR+\n') +
  #labs(title = 'Total Variants\n')+
  #labs(title = 'C-T/G-A Variants\n')+
  #labs(title = 'C-G/G-C Variants\n')+
  #labs(title = 'C-A/G-T Variants\n')+
  #labs(title = 'T-C/A-G Variants\n')+
  #labs(title = 'T-G/A-C Variants\n')+
  #labs(title = 'T-A/A-T Variants\n')+
# s12 
  #xlab('\nVAF Individual 2') + ylab('VAF Mean\n') +
  #labs(title = 'Total Variants\n')+
  #labs(title = 'C-T/G-A Variants\n')+
  #labs(title = 'C-G/G-C Variants\n')+
  #labs(title = 'C-A/G-T Variants\n')+
  #labs(title = 'T-C/A-G Variants\n')+
  #labs(title = 'T-G/A-C Variants\n')+
  #labs(title = 'T-A/A-T Variants\n')+
# s13 
  #xlab('\nVAF Individual 19') + ylab('VAF Mean\n') +
  #labs(title = 'Total Variants\n')+
  #labs(title = 'C-T/G-A Variants\n')+
  #labs(title = 'C-G/G-C Variants\n')+
  #labs(title = 'C-A/G-T Variants\n')+
  #labs(title = 'T-C/A-G Variants\n')+
  #labs(title = 'T-G/A-C Variants\n')+
  #labs(title = 'T-A/A-T Variants\n')+

# 9.20.2017 Analysis
  #xlab('\nVAF Individual 28') + ylab('VAF Mean\n') +

# 9.21.2017 Analysis
  #xlab('\nVAF Individual 19') + ylab('VAF Mean\n') +
  #labs(title = 'Muliplier = 1x\n')+

# 10.25.2017 Analysis
  #xlab('\nVAF Individual 2 Experiment 2') + ylab('VAF Mean\n') +
  #labs(title = 'T-N/A-N Variants')+

# 1.22.2018 Down Syndrome Analysis
  #xlab('\nMean VAF T21') + ylab('Mean VAF D21\n') +
  xlab('\n71.fastq 220102MC') + ylab('72.fastq 220102RUL\n') +
  labs(title = 'Lung Brushing Samples')+

  #geom_smooth(method=lm, se=TRUE, size=3, colour='red')+ # regression line
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

