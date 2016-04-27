taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1] 
ntasks <- as.numeric(commandArgs(trailingOnly=TRUE))[2] 
tname <- commandArgs(trailingOnly=TRUE)[3]
targetpt <- commandArgs(trailingOnly=TRUE)[4]
DST <- commandArgs(trailingOnly=TRUE)[5]

location<-"../scratch/"
tag <- "20160313p"

rDSTall <- ifelse(targetpt=="DS", TRUE, FALSE)
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
{alldrout <- rbind(alldrout, read.csv(paste0(location,"DRcalibration_",currenttag,".",i,".csv"))); i <- i+1} #saved results from dr sampling runs at time 0

tolerance <- 1.5
drout <- alldrout[alldrout[,"rrinc"]/alldrout[,"inc"] > 1/tolerance*alldrout[,"targetdr"] & alldrout[,"rrinc"]/alldrout[,"inc"] < tolerance*alldrout[,"targetdr"], ]  #within 3fold if rr incident fraction target

ilimits <- ceiling(seq(0,nrow(drout), length=ntasks+1))
drout <- drout[(ilimits[taskid]+1):ilimits[taskid+1], ]

currenttag <- paste0(currenttag,".",taskid)

r <- (nrow(read.csv(paste0(location,"TRPwideoutput_", targetpt,DST,"_",currenttag,".csv")))+1):(nrow(drout))

evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, rows=r, targetpt=targetpt, DST=DST, tag=currenttag, rDSTall=rDSTall, location=location, HIV="nonHIV") # can also specify ids and idr to run just a subset of drout
