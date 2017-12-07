# cool little script that checks if ggplot3 is already installed and if not, installs it
list.of.packages <- c('ggplot2')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dep=TRUE)
library(ggplot2)

# Laptop
setwd('/media/alex/Extra/Dropbox/Code/FERMI/baseChangeBias')
# Desktop
#setwd('/home/alex/Dropbox/Code/FERMI/baseChangeBias')
oncoGenes <- read.table('oncoGenes.txt', header = TRUE)
seq=c("A", "C", "T", "G")
df=data.frame('Ref'=rep(seq, each=4), 'Var'=rep(seq, 4), 'Obs'=oncoGenes$Obs)
ggplot(data=df) + aes(x=Var, y=Obs, title='OncoGenes') +
  xlab('Variant')+ ylab('# Observations')+
  geom_bar(stat='identity', position="dodge", fill="grey", color="black")  +
  facet_grid(~Ref)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(size = 50, colour="black"))+
  theme(axis.text.y = element_text(size = 50, colour="black"))+
  theme(axis.ticks = element_line(colour = "black", size = 2))+
  theme(plot.title = element_text(size = 50))+
  theme(plot.title = element_text(hjust = 0.5))+ # center title
  theme(axis.title = element_text(size = 50))+ # change label size 
  theme(strip.text = element_text(size=50)) # change header size

oncoSites <- read.table('oncoSites.txt', header = TRUE)
seq=c("A", "C", "T", "G")
df=data.frame('Ref'=rep(seq, each=4), 'Var'=rep(seq, 4), 'Obs'=oncoSites$Obs)
ggplot(data=df) + aes(x=Var, y=Obs, title='OncoSites') + 
  xlab('Variant')+ ylab('# Observations')+
  geom_bar(stat='identity', position="dodge", fill="grey", color="black")  +
  facet_grid(~Ref)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(size = 50, colour="black"))+
  theme(axis.text.y = element_text(size = 50, colour="black"))+
  theme(axis.ticks = element_line(colour = "black", size = 2))+
  theme(plot.title = element_text(size = 50))+
  theme(plot.title = element_text(hjust = 0.5))+ # center title
  theme(axis.title = element_text(size = 50))+ # change label size 
  theme(strip.text = element_text(size=50)) # change header size

TIIIRegions <- read.table('TIIIRegions.txt', header = TRUE)
seq=c("A", "C", "T", "G")
df=data.frame('Ref'=rep(seq, each=4), 'Var'=rep(seq, 4), 'Obs'=TIIIRegions$Obs)
ggplot(data=df) + aes(x=Var, y=Obs, title='TIIIRegions') + 
  xlab('Variant')+ ylab('# Observations')+
  geom_bar(stat='identity', position="dodge", fill="grey", color="black")  +
  facet_grid(~Ref)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(size = 50, colour="black"))+
  theme(axis.text.y = element_text(size = 50, colour="black"))+
  theme(axis.ticks = element_line(colour = "black", size = 2))+
  theme(plot.title = element_text(size = 50))+
  theme(plot.title = element_text(hjust = 0.5))+ # center title
  theme(axis.title = element_text(size = 50))+ # change label size 
  theme(strip.text = element_text(size=50)) # change header size
