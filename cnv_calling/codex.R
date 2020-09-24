library(CODEX)
dirPath <- "/home/gordeeva/./comparasion_study/exome_data/"
bamFile <- list.files(dirPath, pattern ='*.bam$')
bamdir <- file.path(dirPath, bamFile)
sampname <- as.matrix(substr(bamFile,1,7))
bedFile <- "/home/gordeeva/./comparasion_study/20130108.exome.targets.bed"

finalcall<-matrix(nrow=0,ncol=13,dimnames =list(c(),
                                                c("sample_name", "chr", "cnv","st_bp", "ed_bp",
                                                  "length_kb", "st_exon", "ed_exon"," raw_cov",
                                                  "norm_cov","copy_no","lratio","mBIC")))
for(chr in 1:22){

  bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile, 
                       sampname = sampname, projectname = "CODEX",chr)
  bamdir <- bambedObj$bamdir; sampname <- bambedObj$sampname
  ref <- bambedObj$ref; projectname <- bambedObj$projectname; chr <- bambedObj$chr

  coverageObj <- getcoverage(bambedObj, mapqthres = 20)
  Y <- coverageObj$Y; readlength <- coverageObj$readlength

  gc <- getgc(chr, ref)
  mapp <- getmapp(chr, ref)
  qcObj <- qc(Y, sampname, chr, ref, mapp, gc, cov_thresh = c(20, 4000),
            length_thresh = c(20, 2000), mapp_thresh = 0.9, gc_thresh = c(20, 80))
  Y_qc <- qcObj$Y_qc; sampname_qc <- qcObj$sampname_qc; gc_qc <- qcObj$gc_qc
  mapp_qc <- qcObj$mapp_qc; ref_qc <- qcObj$ref_qc; qcmat <- qcObj$qcmat
  normObj <- normalize(Y_qc, gc_qc, K = 1:5)
  Yhat <- normObj$Yhat; AIC <- normObj$AIC; BIC <- normObj$BIC
  RSS <- normObj$RSS; K <- normObj$K
  optK = K[which.max(BIC)]
  finalcall <-rbind(finalcall, segment(Y_qc, Yhat, optK = optK, K = K, sampname_qc,
                     ref_qc, chr, lmax = 200, mode = "integer"))
}


write.table(finalcall, file = paste('/home/gordeeva/./comparasion_study/calling_tools/codex/CODEX_frac.txt', sep=''), sep='\t', quote=FALSE, row.names=FALSE)