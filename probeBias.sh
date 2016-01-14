'''
This script looks at final .vcf files and for a given probe prints the coverage.
This allows probe bias to be seen as it will show if certain probes are being
represented more or not
'''

INPUT=/home/alex/Desktop/12.22.2015_HiSeq/vcfs/*vcf
OUTPUT=/home/alex/Desktop/12.22.2015_HiSeq/output.txt

declare -a arr=("5073770" "7577539" "7578402" "7577119" "115256529" "115258747" "115258744" "534287" "534288" "534289" "25398284" "25380275" "25457242" "25457243" "209113112" "209113113" "90631934" "90631838" "48649700" "198266832" "198266833" "198266834" "198266835")
for i in "${arr[@]}"
do
    echo "$i"

    awk -F'[;\t]' -v var="$i" '$2 ~ var{sub(/DP=/,"",$15); print $15}' $INPUT >> "$i".txt
done
