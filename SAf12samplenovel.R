taskid <- 2
tname <- "SouthAfrica"
targetpt <- "DR"
DST <- "DSTall"
# pessimistic<-commandArgs(trailingOnly=TRUE)[6]
location<-""

rDSTall <- ifelse(targetpt=="DS", TRUE, FALSE)


tag <- "20160111"

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

drout <- alldrout[alldrout$idr == taskid,]
tolerance <- 1.5
drout <- drout[drout[,"rrinc"]/drout[,"inc"] > 1/tolerance*drout[,"targetdr"] & drout[,"rrinc"]/drout[,"inc"] < tolerance*drout[,"targetdr"], ]  #within 3fold if rr incident fraction target

evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, targetpt=targetpt, DST=DST, tag=currenttag, rDSTall=rDSTall, location=location) # can also specify ids and idr to run just a subset of drout
