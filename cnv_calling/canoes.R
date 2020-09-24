#bedtools multicov -bams `cat /home/gordeeva/./comparasion_study/bam_list` -bed /home/gordeeva/./comparasion_study/20130108.exome.targets.bed -q 20 > /home/gordeeva/./comparasion_study/calling_tools/canoes/canoes.reads.txt

#java -jar /home/gordeeva/tools/GATK/3.5-0-g36282e4/GenomeAnalysisTK.jar -T GCContentByInterval -L /home/gordeeva/./comparasion_study/20130108.exome.targets.bed -R /home/gordeeva/human_genome/hs37d5.fa -o /home/gordeeva/./comparasion_study/calling_tools/canoes/gc.txt

gc <- read.table("/home/gordeeva/./comparasion_study/calling_tools/canoes/gc.txt")$V2
canoes.reads <- read.table("/home/gordeeva/./comparasion_study/calling_tools/canoes/canoes.reads.txt")

sample.names <- c('NA06986.bam','NA06989.bam','NA07051.bam','NA07347.bam','NA10851.bam',
                  'NA11843.bam','NA12340.bam','NA12761.bam','NA12878.bam','NA18959.bam',
                  'NA18960.bam','NA18999.bam')
names(canoes.reads) <- c("chromosome", "start", "end","exon", sample.names)

target <- seq(1, nrow(canoes.reads))
canoes.reads <- cbind(target, gc, canoes.reads)

source("/home/gordeeva/tools/canoes/CANOES.R")

xcnv.list <- vector('list', length(sample.names)) 
for (i in 1:length(sample.names)){
  xcnv.list[[i]] <- CallCNVs(sample.names[i], canoes.reads)
}
xcnvs <- do.call('rbind', xcnv.list)

write.table(xcnvs,paste0("/home/gordeeva/./comparasion_study/calling_tools/canoes/results.csv"),col.names = T,row.names = F,quote = F)