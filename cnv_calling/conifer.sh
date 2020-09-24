#!/bin/bash

cd /home/gordeeva/./comparasion_study/calling_tools/conifer
grep -v "Y" /home/gordeeva/./comparasion_study/20130108.exome.targets.bed > probes.txt 

for sample in NA06986 NA06989 NA07051 NA07347 NA11843 NA12340 NA12761 NA12878 NA18959 NA18960 NA18999
do
python /home/gordeeva/tools/conifer_v0.2.2/conifer.py  rpkm --probes probes.txt --input /home/gordeeva/./comparasion_study/exome_data/$sample.bam --output RPKM/$sample.rpkm.txt
done

python /home/gordeeva/tools/conifer_v0.2.2/conifer.py analyze --probes probes.txt --rpkm_dir ./RPKM/ --output analysis.hdf5 --svd 2 --write_svals singular_values.txt --plot_scree screeplot.png --write_sd sd_values.txt

python/home/gordeeva/tools/conifer_v0.2.2/conifer.py call \
  --input analysis.hdf5 \
  --output calls.txt    \
   --threshold 0.95