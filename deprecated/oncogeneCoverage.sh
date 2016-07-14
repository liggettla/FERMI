#The purpose of this script is to look at the gene coverage of the oncogene
#sites to look for any coverage differences in samples. This script outputs
#Dp values to file that can then be processed by oncogeneCovverage.py

grep --color=always '50737*\|75775*\|75771*\|1152565*\|1152587*\|1152587*\|534287\|5342*\|5342*\|253982*\|253802*\|1061972*\|1061972*\|1061972*\|1061972*\|1061551*\|1061551*\|1061551*\|254572*\|254572*\|2091131*\|2091131*\|906319*\|906318\|486497*' 5_finalOutput.vcf | less | cut -f 8 | awk -F \; '{sub(/DP=/,"")} {print $8}' >> coverageOutput.txt
