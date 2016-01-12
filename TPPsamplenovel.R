taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1] #should have same array as dr runs
tname <- commandArgs(trailingOnly=TRUE)[2]
targetpt <- commandArgs(trailingOnly=TRUE)[3]
DST <- commandArgs(trailingOnly=TRUE)[4]
rDSTall<-commandArgs(trailingOnly=TRUE)[5]
location<-"../scratch/"


tag <- "20160111"
currenttag <- paste0(tname,"_",tag,".",taskid)
if(rDSTall==TRUE) currenttag <- paste0("rDSTall.",currenttag)

source("TPPmat.R")

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
values <- set.values()
mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- set.novelvalues()$elementnames

drout <- screendrout(drout_filename=paste0(location,"DRcalibration_",currenttag,".csv"), tolerance=1.5) # limits output to those simulations with DR fraction of incidence in targetepi range

startrow <- max(read.csv(paste0(location,"TRPwideoutput_",targetpt,DST,"_",currenttag,".csv"))$inew) +1

evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, rows=startrow:nrow(drout), targetpt=targetpt, DST=DST, tag=currenttag, rDSTall=rDSTall, location=location) # can also specify ids and idr to run just a subset of drout
