'''
This script looks at final .vcf files and for a given probe prints the coverage.
This allows probe bias to be seen as it will show if certain probes are being
represented more or not
'''

INPUT=/home/alex/Desktop/12.22.2015_HiSeq/vcfs/*vcf
OUTPUT=/home/alex/Desktop/12.22.2015_HiSeq/output.txt

declare -a arr=("50737" "75775" "75784" "75771" "1152565" "1152587" "1152587" "53428" "53428" "53428" "253982" "253802" "254572" "254572" "2091131" "2091131" "906319" "906318" "486497" "1982668" "1982668" "1982668" "1982668")
for i in "${arr[@]}"
do
    echo "$i"

    awk -F'[;\t]' -v var="$i" '$2 ~ var{sub(/DP=/,"",$15); print $15}' $INPUT >> "$i".txt
done
