#!/bin/bash


###config file for ExonDel software
#reference .bed file
#bedfile=/home/gordeeva/./comparasion_study/calling_tools/exomedepth/refseq.bed ##
#reference .gtf file
#refseq=/home/gordeeva/./comparasion_study/calling_tools/exomedepth/exondel/refseq.gtf
#reference .fa file
#reffa=/home/gordeeva/human_genome/hs37d5.fa


cd /home/gordeeva/./comparasion_study/calling_tools/exomedepth/

perl /home/gordeeva/tools/ExonDel-master/ExonDel.pl -i bam_list -c /home/gordeeva/tools/ExonDel-master/ExonDel.cfg \
    -o ./results