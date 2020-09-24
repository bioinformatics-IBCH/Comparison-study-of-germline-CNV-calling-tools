#!/bin/bash
cd /home/gordeeva/./comparasion_study/exome_data/
python /home/gordeeva/tools/CONTRA.v2.0.8/baseline.py  --target /home/gordeeva/./comparasion_study/20130108.exome.targets.bed --files NA06986.bam NA06989.bam NA07051.bam NA07347.bam NA11843.bam NA12340.bam NA12761.bam NA18959.bam NA18960.bam NA18999.bam --output /home/gordeeva/./comparasion_study/calling_tools/contra/baseline/ --name baseline_test


python /home/gordeeva/tools/CONTRA.v2.0.8/contra.py --target /home/gordeeva/./comparasion_study/20130108.exome.targets.bed  --test NA12878.bam -c  /home/gordeeva/./comparasion_study/calling_tools/contra/baseline/baseline_test.pooled2_TRIM0.2.txt  --bed -o /home/gordeeva/./comparasion_study/calling_tools/contra/results