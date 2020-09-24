library(exomeCopy)
target.file <- "/home/gordeeva/./comparasion_study/20130108.exome.targets.bed"
bam.files <- c("/home/gordeeva/./comparasion_study/exome_data/NA06986.bam",
               "/home/gordeeva/./comparasion_study/exome_data/NA06989.bam",
               "/home/gordeeva/./comparasion_study/exome_data/NA07051.bam",
               "/home/gordeeva/./comparasion_study/exome_data/NA07347.bam",
               "/home/gordeeva/./comparasion_study/exome_data/NA11843.bam",
               "/home/gordeeva/./comparasion_study/exome_data/NA12340.bam",
               "/home/gordeeva/./comparasion_study/exome_data/NA12761.bam",
               "/home/gordeeva/./comparasion_study/exome_data/NA12878.bam",
               "/home/gordeeva/./comparasion_study/exome_data/NA18959.bam",
               "/home/gordeeva/./comparasion_study/exome_data/NA18960.bam",
               "/home/gordeeva/./comparasion_study/exome_data/NA18999.bam")
sample.names <- c("NA06986","NA06989","NA07051","NA07347","NA11843","NA12340",
                  "NA12761","NA12878","NA18959","NA18960","NA18999")

reference.file <- "/home/gordeeva/human_genome/hs37d5.fa"


target.df <- read.delim(target.file, header = FALSE)
target <- GRanges(seqname = target.df[, 1], IRanges(start = target.df[,2] + 1, 
                                                    end = target.df[, 3]))
counts <- RangedData(space = seqnames(target), ranges = ranges(target))
for (i in 1:length(bam.files)) {
  counts[[sample.names[i]]] <- countBamInGRanges(bam.files[i],target)}
counts[["GC"]] <- getGCcontent(target, reference.file)
counts[["GC.sq"]] <- counts$GC^2
counts[["bg"]] <- generateBackground(sample.names, counts, median)
counts[["log.bg"]] <- log(counts$bg + 0.1)
counts[["width"]] <- width(counts)


chr<-c("1","10","11","12","13","14","15","16","17","18","19",
       "2","20","21","22","3","4","5","6","7","8","9")


fit.list <-   lapply(chr, function(seq.name) {
    exomeCopy(counts[seq.name],
              sample.name=sample.names[8], X.names = c("log.bg", "GC","GC.sq", "width"),
              S = 0:4, d = 2)
    })



compiled.segments<-c(copyCountSegments(fit.list[[1]]),copyCountSegments(fit.list[[2]]),
                     copyCountSegments(fit.list[[3]]),copyCountSegments(fit.list[[4]]),
                     copyCountSegments(fit.list[[5]]),copyCountSegments(fit.list[[5]]),
                     copyCountSegments(fit.list[[7]]),copyCountSegments(fit.list[[8]]),
                     copyCountSegments(fit.list[[9]]),copyCountSegments(fit.list[[10]]),
                     copyCountSegments(fit.list[[11]]),copyCountSegments(fit.list[[12]]),
                     copyCountSegments(fit.list[[13]]),copyCountSegments(fit.list[[14]]),
                     copyCountSegments(fit.list[[15]]),copyCountSegments(fit.list[[16]]),
                     copyCountSegments(fit.list[[17]]),copyCountSegments(fit.list[[18]]),
                     copyCountSegments(fit.list[[19]]),copyCountSegments(fit.list[[20]]),
                     copyCountSegments(fit.list[[21]]),copyCountSegments(fit.list[[22]]))
CNV.segments <- compiled.segments[compiled.segments$copy.count != 2, ]
CNV.segments <- CNV.segments[CNV.segments$nranges > 5,]

CNV.segments<-as(CNV.segments, "data.frame")
write.table(CNV.segments,paste0("/home/gordeeva/./comparasion_study/calling_tools/exomecopy/results"),col.names = T,row.names = F,quote = F,sep="\t")
