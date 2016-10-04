library(ggplot2)

setwd('/home/alex/Dropbox/Degregori_Lab/Experiments/10.4.2016/A1-R2.fastq')
oncoGenes <- read.table('oncoGenes.txt', header = TRUE)
seq=c("A", "C", "T", "G")
df=data.frame('Ref'=rep(seq, each=4), 'Var'=rep(seq, 4), 'Obs'=oncoGenes$Obs)
ggplot(data=df) + aes(x=Var, y=Obs, title='A1-R2 OncoGenes') +
  geom_bar(stat='identity', position="dodge", fill="grey", color="black")  +
  facet_grid(~Ref)

oncoSites <- read.table('oncoSites.txt', header = TRUE)
seq=c("A", "C", "T", "G")
df=data.frame('Ref'=rep(seq, each=4), 'Var'=rep(seq, 4), 'Obs'=oncoSites$Obs)
ggplot(data=df) + aes(x=Var, y=Obs, title='A1-R2 OncoSites') + 
  geom_bar(stat='identity', position="dodge", fill="grey", color="black")  +
  facet_grid(~Ref)

TIIIRegions <- read.table('TIIIRegions.txt', header = TRUE)
seq=c("A", "C", "T", "G")
df=data.frame('Ref'=rep(seq, each=4), 'Var'=rep(seq, 4), 'Obs'=TIIIRegions$Obs)
ggplot(data=df) + aes(x=Var, y=Obs, title='A1-R2 TIIIRegions') + 
  geom_bar(stat='identity', position="dodge", fill="grey", color="black")  +
  facet_grid(~Ref)