#!/bin/bash

cd /home/gordeeva/./comparasion_study/calling_tools/fishingcnv/

java -Xmx6g -Xms2g -jar GenomeAnalysisTK.jar \
    -T DepthOfCoverage -I /home/gordeeva/./comparasion_study/exome_data/NA#.bam \
    -L /home/gordeeva/./comparasion_study/20130108.exome.targets.bed \
    -R /home/gordeeva/human_genome/hs37d5.fa\
    -o output NA#.coverage \
    --minMappingQuality 15 --minBaseQuality 10 --omitDepthOutputAtEachBase \
    --logging_level INFO --summaryCoverageThreshold 5 --summaryCoverageThreshold 7 \
    --summaryCoverageThreshold 10 --summaryCoverageThreshold 15 \
    --summaryCoverageThreshold 20 --summaryCoverageThreshold 30 --summaryCoverageThreshold 50

java -jar /home/gordeeva/tools/FishingCNV_2.1_pipeline/FishingCNV.jar -cc -c NA#.coverage.sample_interval_summary \
-b /home/gordeeva/./comparasion_study/20130108.exome.targets.bed -o NA#.rpkm

java -jar /home/gordeeva/tools/FishingCNV_2.1_pipeline/FishingCNV.jar -p -rpkm NA06986.rpkm NA06989.rpkm NA07051.rpkm NA07347.rpkm \ 
NA11843.rpkm NA12340.rpkm NA12761.rpkm NA18959.rpkm NA18960.rpkm NA18999.rpkm -o controls.ctr

Rscript /home/gordeeva/tools/FishingCNV_2.1_pipeline/FishingCNV.jar -c controls.ctr.complete -v -s NA12878.rpkm
 -o results -pca