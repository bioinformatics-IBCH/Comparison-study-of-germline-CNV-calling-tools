#!/bin/bash

cd /home/gordeeva/tools/EXCAVATOR2_Package_v1.1.2

#cat SourceTarget.txt 
#data/ucsc.hg19.bw /home/gordeeva/human_genome/hs37d5.fa

perl TargetPerla.pl  SourceTarget.txt /home/gordeeva/./comparasion_study/20130108.exome.targets.bed  1KG 50000 hg19

#PrepareFile.txt Tab-delimited)
#/home/gordeeva/./comparasion_study/exome_data/NA06986.bam    /home/gordeeva/./comparasion_study/calling_tools/excavator/DataPrepare/NA0698    NA06986

perl EXCAVATORDataPrepare.pl PrepareFile.txt --processors 6 --target 1KG --assembly hg19

#AnalysisFile.txt Tab-delimited)
#C1 /home/gordeeva/./comparasion_study/calling_tools/excavator/DataPrepare/NA06986    NA06986
#C2 /home/gordeeva/./comparasion_study/calling_tools/excavator/DataPrepare/NA06989    NA06989
#...
#T1 /home/gordeeva/./comparasion_study/calling_tools/excavator/DataPrepare/NA12878    NA12878

perl EXCAVATORDataAnalysis.pl AnalysisFile.txt --processors 1 --target 1KG --assembly hg19  --output /home/gordeeva/./comparasion_study/calling_tools/excavator/results --mode pooling