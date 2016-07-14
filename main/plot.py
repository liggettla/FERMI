# This does nothing at the moment but will be use to automatically generate plots

'''
Nothing below this has yet been implemented
#######################
#Output Var Freq Table#
#######################
outputDir = outputDir + '/'
vcfFile = outputDir + '/finalOutput.vcf'
system('cp %s ./' % (vcfFile))
system('bash newIdentifyVars.sh')
system('mv finalOutput.vcf %s' (outputDir))

###################
#Plot Allele Freqs#
###################
system('Rscript plotVarFreq.R')
system('mv allelefreqs.txt %s' (outputDir))
system('mv allelefreqs.jpg %s' (outputDir))


#this should output a plot automatically, but is not runtime tested
#does not output log10() of plot yet
#it may be best to just system(Rscript plotting.R) somehow
if False:
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
r = ro.r
r.setwd(outputDir)
f = r('read.table("allelefreqs.txt", header = FALSE)')
grdevices = importr('grDevices')
grdevices.png(file="alleleFreq.png", width=800, height=500)
r.hist(f[0], breaks=100, main = '5 Reads', xlab='Variant Freq', ylab='# Vars')
grdevices.dev_off()
'''
