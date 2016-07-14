#The purpose of this script is to count the number of Onc/TIII variants in a sample

#These are the probed oncogenic bp locations
grep '5073770\|7577539\|7577119\|115256529\|115258747\|115258744\|534287\|534288\|534289\|25398284\|25380275\|106197266\|106197267\|106197268\|106197269\|106155172\|106155173\|106155174\|25457242\|25457243\|209113112\|209113113\|90631934\|90631838\|48649700' *vcf | wc -l

#These are all oncogene locations but the above grep will get just the locs
#that are oncogenic and can be subtracted from the following grep in order
#to get just those locations that are non-oncogenic
grep --color=always '50737*\|75775*\|75771*\|1152565*\|1152587*\|1152587*\|534287\|5342*\|5342*\|253982*\|253802*\|1061972*\|1061972*\|1061972*\|1061972*\|1061551*\|1061551*\|1061551*\|254572*\|254572*\|2091131*\|2091131*\|906319*\|906318\|486497*' *vcf | wc -l
