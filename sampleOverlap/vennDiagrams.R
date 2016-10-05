# cool little script that checks if ggplot2 is already installed and if not, installs it
list.of.packages <- c('VennDiagram')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dep=TRUE)
library(VennDiagram)
setwd("~/Desktop")

jpeg('C1-D1.jpg')
# venn(sample1, sample2, overlap, ....)
draw.pairwise.venn(3099, 3040, 2629, category = c("C1 F34", "D1 F34"), lty = rep("blank", 
  2), fill = c("blue", "green"), alpha = rep(0.5, 2), cat.pos = c(193, 
  150), cat.dist = rep(0.05, 2), scaled = TRUE)
dev.off()

jpeg('E1-F1.jpg')
# venn(sample1, sample2, overlap, ....)
draw.pairwise.venn(3134, 2987, 2599, category = c("E1 F41", "F1 F41"), lty = rep("blank", 
  2), fill = c("blue", "green"), alpha = rep(0.5, 2), cat.pos = c(193, 
  150), cat.dist = rep(0.05, 2), scaled = TRUE)
dev.off()

jpeg('G1-H1.jpg')
# venn(sample1, sample2, overlap, ....)
draw.pairwise.venn(3092, 3217, 2665, category = c("G1 F34", "H1 F34"), lty = rep("blank", 
  2), fill = c("blue", "green"), alpha = rep(0.5, 2), cat.pos = c(193, 
  150), cat.dist = rep(0.05, 2), scaled = TRUE)
dev.off()

