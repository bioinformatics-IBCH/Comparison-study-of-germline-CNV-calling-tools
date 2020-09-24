#!/bin/bash

/home/gordeeva/tools/cnvkit/cnvkit.py batch /home/gordeeva/./comparasion_study/exome_data/NA12878.bam \
    --normal /home/gordeeva/./comparasion_study/exome_data/NA06986.bam \
             /home/gordeeva/./comparasion_study/exome_data/NA06989.bam \
             /home/gordeeva/./comparasion_study/exome_data/NA07051.bam \ 
             /home/gordeeva/./comparasion_study/exome_data/NA07347.bam \
             /home/gordeeva/./comparasion_study/exome_data/NA11843.bam \
             /home/gordeeva/./comparasion_study/exome_data/NA12340.bam \
             /home/gordeeva/./comparasion_study/exome_data/NA12761.bam \
             /home/gordeeva/./comparasion_study/exome_data/NA12878.bam \
             /home/gordeeva/./comparasion_study/exome_data/NA18959.bam \
             /home/gordeeva/./comparasion_study/exome_data/NA18960.bam \
             /home/gordeeva/./comparasion_study/exome_data/NA18999.bam \
    --targets /home/gordeeva/./comparasion_study/20130108.exome.targets.bed \
    --annotate /home/gordeeva/tools/cnvkit/refFlat.txt \
    --fasta /home/gordeeva/human_genome/hs37d5.fa \
    --access /home/gordeeva/tools/cnvkit/data/access-5kb-mappable.hg19.bed \
    --output-reference /home/gordeeva/tools/cnvkit/my_reference.cnn\
    --output-dir /home/gordeeva/./comparasion_study/calling_tools/cnvkit/ \
    --diagram --scatter