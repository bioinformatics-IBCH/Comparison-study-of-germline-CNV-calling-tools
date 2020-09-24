library(GenomicFeatures)
library(DBI)

txdb <- makeTxDbFromGFF(file="/home/gordeeva/human_genome/gencode.v19.annotation.gff3")

exons.list.per.gene <- exonsBy(txdb,by="gene")
k<- lapply(exons.list.per.gene, function(x) reduce(x))
k<-unlist(GRangesList(k))
values(k) <- DataFrame(gene_id = names(k))
exons_df = as.data.frame(k,row.names = seq(1, 324422))
exons_df <-exons_df[order(exons_df$seqnames, exons_df$start),]
exons_df$seqnames<-str_split_fixed(exons_df$seqnames, "chr", 2)[,2] 

write.csv(exons_df, file = "/home/gordeeva/human_genome/gencode_exons", row.names=FALSE)