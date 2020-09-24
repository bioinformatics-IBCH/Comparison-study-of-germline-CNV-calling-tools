library(Biostrings)
library(IRanges)
library(ExomeDepth)

bed <- read.table("/home/gordeeva/./comparasion_study/20130108.exome.targets.bed",sep="\t",header = FALSE) #path of BED-File exons 
#exome data- file with paths to processed BAM files
samples.bam <- scan("/home/gordeeva/./comparasion_study/exome_data", what="", sep="\n")


fasta <- "/home/gordeeva/human_genome/hs37d5.fa"

read.count <- getBamCounts(bed.frame = bed,
                          bam.files = samples.bam,
                          include.chr = FALSE,
                          referenceFasta = fasta)

save(reads.counts, file="/home/gordeeva/./comparasion_study/calling-tools/exomedepth/counts.RData")

rc <- as.data.frame(read.count)

my.test<-rc$NA12878.bam
my.reference.set <-as.matrix(rc[,6:ncol(norma),select=-NA12878.bam])
my.choice<-select.reference.set(test.counts= my.test,                                  
                                reference.counts= my.reference.set, 
                                bin.length = (rc$end-rc$start)/1000,
                                n.bins.reduced=10000)
write.table(head(my.choice$summary.stats,10),"/home/gordeeva/./comparasion_study/reference.txt",col.names = T,row.names = F,quote = F)