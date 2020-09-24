#!/bin/bash

/home/gordeeva/tools/DeAnnCNV/PreprocessFiles/run.sh bam_list \
    /home/gordeeva/human_genome/hs37d5.fa \
    /home/gordeeva/./comparasion_study/20130108.exome.targets.bed \
    /home/gordeeva/./comparasion_study/calling_tools/deanncnv
    
# load tar.gz to https://mcg.ustc.edu.cn/bsc/cnv/upload.html