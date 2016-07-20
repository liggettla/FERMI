library(VennDiagram)

# Sample repeatability between Cord samples 35 and 36 from 3.21.2016 experiment
# represents 84% and 62% overlap
#grid.newpage()
jpeg('rplot.jpg')
draw.pairwise.venn(3063, 4158, 2572, category = c("Cord 1 (Sample 35)", "Cord 2 (Sample 36)"), lty = rep("blank", 
  2), fill = c("blue", "green"), alpha = rep(0.5, 2), cat.pos = c(193, 
  150), cat.dist = rep(0.05, 2), scaled = TRUE)
dev.off()

# Sample repeatability between 293 samples 33 and 34 from 3.21.2016 experiment
grid.newpage()
draw.pairwise.venn(3063, 4158, 2572, category = c("Cord 1 (Sample 35)", "Cord 2 (Sample 36)"), lty = rep("blank", 
  2), fill = c("blue", "green"), alpha = rep(0.5, 2), cat.pos = c(193, 
  150), cat.dist = rep(0.05, 2), scaled = TRUE)

# Sample repeatability between B Cells (sample 31) and Myeloid Cells (sample 32) from 3.21.2016 experiment
# represents a 67% and 73% overlap
grid.newpage()
draw.pairwise.venn(2496, 2321, 1683, category = c("B Cells (Sample 31)", "Myeloid Cells (Sample 32)"), lty = rep("blank", 
  2), fill = c("blue", "green"), alpha = rep(0.5, 2), cat.pos = c(-20, 
  15), cat.dist = rep(0.05, 2), scaled = TRUE)
