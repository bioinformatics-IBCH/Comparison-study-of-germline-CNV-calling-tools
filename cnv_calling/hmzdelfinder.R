library(RCurl)
library(data.table)
library(gdata)
library(parallel)
library(Hmisc)
library(matrixStats)
library(DNAcopy)
library(GenomicRanges)
library(Rsubread) 

eval( expr = parse( text = getURL("https://raw.githubusercontent.com/BCM-Lupskilab/HMZDelFinder/master/src/HMZDelFinder.R") ))
?mainDir <- '/home/gordeeva/./comparasion_study/calling_tools/hmzdelfinder'


# set/create other paths and identifiers
bedFile <- "/home/gordeeva/./comparasion_study/20130108.exome.targets.bed" 
pathToBams <- "/home/gordeeva/./comparasion_study/exome_data/" 
bamFiles <- paste0(pathToBams, dir(pathToBams, "bam$"))
rpkmDir <-  paste0(mainDir, "rpkm/" , sep=""); if (!file.exists(plotsDir)){dir.create(plotsDir)}  
sampleNames <- sapply(strsplit(dir(pathToBams, "bam$"), "[/\\.]"), function(x){x[length(x)-1]})
calcRPKMsFromBAMs(bedFile,bamFiles , sampleNames, rpkmDir,4)


rpkmFiles <- dir(rpkmDir, "rpkm.txt$")
rpkmFids <- gsub(".rpkm.txt", "", rpkmFiles) 
rpkmPaths <- paste0(rpkmDir, rpkmFiles)
aohDir <- paste0(mainDir, "AOH/" , sep=""); if (!file.exists(aohDir)){dir.create(aohDir)} 
aohRDataOut <- paste(mainDir, "AOH/extAOH_small.RData", sep="")
outputDir <- paste0(mainDir, "out/" , sep=""); if (!file.exists(outputDir)){dir.create(outputDir)} # create output directory
plotsDir <- paste0(mainDir, "plots/" , sep=""); if (!file.exists(plotsDir)){dir.create(plotsDir)} # create output plots directory


is_cmg <- FALSE
lowRPKMthreshold <- 0.65
maxFrequency <- 0.05
minAOHsize <- 1000
minAOHsig <- 0.45
mc.cores<-4
vR_id<-"VR"
tR_id<-"DP"
filter <- "PASS"

# running HMZDelFinder
results <- runHMZDelFinder (NULL,NULL, 
                            rpkmPaths,rpkmFids,
                            mc.cores,aohRDataOut,
                            bedFile,
                            lowRPKMthreshold, minAOHsize, minAOHsig,
                            is_cmg,vR_id,tR_id,filter)

# saving results in csv files
write.csv(results$filteredCalls, paste0(outputDir,"hmzCalls.csv"), row.names=F )
