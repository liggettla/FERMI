#This script finds all of the TIII variations that are at an AF=0
#and outputs a count

#This pulls the info lines from all matching variants then counts
#only those that have and AF (column 4) that == 0
grep '115227*\|229041*\|110541*\|112997*\|121167*\|123547*\|124428*\|1397*\|2126*\|2390*\|2593*\|11486*\|92527*\|73379*\|82455*\|85949*' 16_finalOutput.vcf | cut -f 8 | awk -F \; '{sub(/AF=/,"")} $4 != 1 && $4 != 0.5 {print $4}' | wc -l
