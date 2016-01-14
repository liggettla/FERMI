'''
This script looks at final .vcf files and for a given probe prints the coverage.
This allows probe bias to be seen as it will show if certain probes are being
represented more or not
'''

INPUT=/home/alex/Desktop/12.22.2015_HiSeq/vcfs/*vcf
echo 'Probe Loc:' > output.txt
echo '50737' >> output.txt
echo '' >> output.txt
echo 'Coverage:' >> output.txt
awk -F'[;\t]' '$2 ~ /50737../{sub(/DP=/,"",$15); print $15}' $INPUT >> output.txt
