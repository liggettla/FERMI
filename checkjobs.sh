# the purpose of this is to check current jobs running on the cluster
while true
    do echo                                                                                 "----------------------------------------------------------------------------------" 
    #ls -lht /vol3/home/liggettl/sequencingData/1.12.2016/results/2016-05-03_3_0.51/30_TGTGACCA_L008_R1_001.fastq
    echo ""
    bjobs | sort
    echo "----------------------------------------------------------------------------------"
    sleep 1
    echo ""
done
