taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1] 
ntasks <- as.numeric(commandArgs(trailingOnly=TRUE))[2] 
tname <- commandArgs(trailingOnly=TRUE)[3]
targetpt <- commandArgs(trailingOnly=TRUE)[4]
DST <- commandArgs(trailingOnly=TRUE)[5]

location<-"../scratch/"
tag <- "20160304p"; Nsims_ds <- 250

rDSTall <- ifelse(targetpt=="DS", TRUE, FALSE)
ilimits <- ceiling(seq(0,Nsims_ds, length=ntasks+1))
currenttag <- paste0(tname,"_",tag)
if(rDSTall==TRUE) currenttag <- paste0("rDSTall.",currenttag)

source("TPPmat.R")

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
values <- set.values()
mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- set.novelvalues()$elementnames

alldrout <- numeric(0)
i <- 1; while(file.exists(paste0(location,"DRcalibration_",currenttag,".",i,".csv")))
{alldrout <- rbind(alldrout, read.csv(paste0(location,"DRcalibration_",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0

currenttag <- paste0(currenttag,".",taskid)

drout <- alldrout[alldrout$ids %in% (ilimits[taskid]+1):ilimits[taskid+1],]
tolerance <- 1.5
drout <- drout[drout[,"rrinc"]/drout[,"inc"] > 1/tolerance*drout[,"targetdr"] & drout[,"rrinc"]/drout[,"inc"] < tolerance*drout[,"targetdr"], ]  #within 3fold if rr incident fraction target
# if (fixavail) drout$fixed.availability <- 0.8

evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, targetpt=targetpt, DST=DST, tag=currenttag, rDSTall=rDSTall, location=location) # can also specify ids and idr to run just a subset of drout
