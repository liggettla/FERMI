# cool little script that checks if ggplot2 is already installed and if not, installs it
list.of.packages <- c("ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dep=TRUE)

library(ggplot2)

# read in the data
vafs <- read.table("outputFile", header = TRUE)
sample1 <- vafs$Sample1
sample2 <- vafs$Sample2

# plot with 95% confidence interval
lm_fit  = lm(sample1 ~ sample2)
x = data.frame(vafs, predict(lm_fit, interval = 'prediction'))

p <- ggplot(x, aes(x=sample1,y=sample2, alpha=0.5)) +
  xlab('Sample 1') +
  ylab('Sample 2') +
  labs(title = 'VAF Repeatability') +
  geom_point() +
  geom_smooth(method = 'lm', aes(fill = 'confidence'), alpha = 0.15) +
  scale_fill_manual('Interval', values = c('green', 'blue'))

jpeg('output.jpg')
print(p)
dev.off()