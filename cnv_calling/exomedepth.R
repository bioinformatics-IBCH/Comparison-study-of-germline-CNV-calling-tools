library(Biostrings)
library(IRanges)
library(ExomeDepth)

read.count <- get(load(paste0("/home/gordeeva/./comparasion_study/calling_tools/exomedepth/counts.RData")))
norma <-as.data.frame(read.count)

my.test <-norma$NA12878.bam

x <- c("NA06986.bam","NA06989.bam","NA07051.bam","NA07347.bam","NA11843.bam","NA12340.bam",
       "NA12761.bam","NA18959.bam","NA18960.bam","NA18999.bam")
my.matrix <- as.matrix( norma[, x])
my.reference.selected <- apply(X = my.matrix, MAR=1,FUN=sum)

all.exons <- new('ExomeDepth', test = my.test, 
                 reference = my.reference.selected,
                 formula ='cbind(test, reference) ~ 1')


all.exons <-CallCNVs(x = all.exons, transition.probability = 10^-4,
                     chromosome =factor(as.character(unlist(control1$space))),
                     start = as.numeric(unlist(control1$start)), 
                     end = as.numeric(unlist(control1$end)),name=paste0("exon",1:nrow(control1)))

data<-as.data.frame(all.exons@CNV.calls)
write.table(data,paste0("/home/gordeeva/./comparasion_study/calling_tools/exomedepth/cnv"),col.names = T,row.names = F,quote = F)

