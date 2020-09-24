#!/bin/bash

#######normalize_coverage
#OUTPUT_DIR=/home/gordeeva/./comparasion_study/calling_tools/clamms
#sample=`echo $1| awk -F"/" '{print $NF}'|awk -F"." '{print $1}'`
#echo $sample

#samtools bedcov -Q 30 /home/gordeeva/./comparasion_study/calling_tools/clamms/windows.bed $1 \
#| awk '{ printf "%s\t%d\t%d\t%.6g\n", $1, $2, $3, $NF/($3-$2); }' \
#> $OUTPUT_DIR/${sample}.coverage.bed

#/home/gordeeva/tools/clamms/normalize_coverage $OUTPUT_DIR/${sample}.coverage.bed \
#$OUTPUT_DIR/windows.bed > $OUTPUT_DIR/${sample}.norm.cov.bed
############

export $CLAMMS_DIR=/home/gordeeva/tools/clamms

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig
/home/gordeeva/tools/bigWigToWig wgEncodeCrgMapabilityAlign75mer.bigWig wgEncodeCrgMapabilityAlign75mer.wig
grep -v '^#' wgEncodeCrgMapabilityAlign75mer.wig | sed 's/^chr//g' > $CLAMMS_DIR/mappability.bed

cd /home/gordeeva/./comparasion_study/calling_tools/clamms
export INSERT_SIZE=250
$CLAMMS_DIR/annotate_windows.sh /home/gordeeva/./comparasion_study/20130108.exome.targets.bed ./home/gordeeva/human_genome/hs37d5.fa $CLAMMS_DIR/mappability.bed $INSERT_SIZE $CLAMMS_DIR/data/clamms_special_regions.hg19.bed > windows.bed

cat bam_list | xargs -P $NUM_PROCESSES --max-args 1 $CLAMMS_DIR/normalize_coverage.sh


###ref.panel.files.txt 
#${sample}.norm.cov.bed     F/M

$CLAMMS_DIR/fit_models ref.panel.files.txt windows.bed > /$CLAMMS_DIR/models.bed
$CLAMMS_DIR/call_cnv NA12878.norm.cov.bed models.bed --sex F > NA12878.cnv.bed