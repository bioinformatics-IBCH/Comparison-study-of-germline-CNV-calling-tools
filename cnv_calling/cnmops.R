source("https://bioconductor.org/biocLite.R")
biocLite("cn.mops")
library(cn.mops)
#https://support.bioconductor.org/p/70037/

BAMFiles<-c("/home/gordeeva/./comparasion_study/exome_data/NA06986.bam",
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

segments <- read.table("/home/gordeeva/./comparasion_study/20130108.exome.targets.bed",
                       sep="\t",as.is=TRUE)
gr<-GRanges(segments[,1],IRanges(segments[,2],segments[,3]))
X<-getSegmentReadCountsFromBAM(BAMFiles,GR=gr)

resCNMOPS<-exomecn.mops(X)
resCNMOPS<-calcIntegerCopyNumbers(resCNMOPS)
CNVs<-as.data.frame(cnvs(resCNMOPS))
write.csv(CNVs,file="/home/gordeeva/./comparasion_study/calling_tools/cnmops/cnmops.csv")