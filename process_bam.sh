#Exome data processing
#In: file with paths to BAM files
#    number of processors


#!/bin/bash
sample_list=$1
n_core=$2
func() {
        file=`echo $1 | awk '{print $0}'`
        sample=`echo $file | awk -F "/" '{print $NF}' |  awk -F "." '{print $1}'`
        echo "Processing "$sample
        #Create directory
        [ -d  exome_data] || mkdir exome_data
        
        echo "Filter low quality reads..."
        samtools view -b -q 11 $file | samtools sort -o exome_data/$sample.sorted.bam

        /home/gordeeva/tools/picard-2.18.2-0/picard MarkDuplicates \
        I= exome_data/$sample.sorted.bam O=exome_data/$sample.bam METRICS_FILE=$sample_output.dup_metrics \
        REMOVE_DUPLICATES=TRUE  VALIDATION_STRINGENCY=SILENT

        samtools index exome_data/$sample.bam && rm exome_data/$sample.sorted.bam

}

export -f func
cat $sample_list | parallel -j"$n_core" -a - func