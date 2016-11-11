#!/usr/bin/python

from parseLine import seqRead


line = 'chr18\t22343095\t.\tG\tA\t13777\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=636;CIGAR=1X;DP=637;DPB=637;DPRA=0;EPP=1384.07;EPPR=5.18177;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=60;NS=1;NUMALT=1;ODDS=875.012;PAIRED=0;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=21725;QR=36;RO=1;RPL=636;RPP=1384.07;RPPR=5.18177;RPR=0;RUN=1;SAF=636;SAP=1384.07;SAR=0;SRF=1;SRP=5.18177;SRR=0;TYPE=snp\tGT:DP:DPR:RO:QR:AO:QA:GL\t1/1:637:637,636:1:36:636:21725:-1951.16,-188.158,0\n'

x = seqRead(line)
print x.chrom()
print x.loc()
print x.wt()
print x.var()
print x.ao()
print x.dp()
print x.af()
