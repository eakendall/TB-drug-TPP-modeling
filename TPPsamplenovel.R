taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1] #should have same array as dr runs
tname <- commandArgs(trailingOnly=TRUE)[2]
targetpt <- commandArgs(trailingOnly=TRUE)[3]
DST <- commandArgs(trailingOnly=TRUE)[4]
#startrow <- as.numeric(commandArgs(trailingOnly=TRUE))[5]

tag <- "20151227"
currenttag <- paste0(tname,"_",tag,".",taskid)

source("TPPmat.R")

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
values <- set.values()
mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- set.novelvalues()$elementnames


drout <- screendrout(drout_filename=paste0("../scratch/DRcalibration_",currenttag,".csv"), tolerance=1.5) # limits output to those simulations with DR fraction of incidence in targetepi range

startrow <- max(read.csv(paste0("../scratch/TRPwideoutput_",targetpt,DST,"_",currenttag,".csv"))$inew) +1

idr <-  1:20 # use for test runs, or if splitting DR tasks into subtasks here

evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout[startrow:nrow(drout), ], idr=idr, targetpt=targetpt, DST=DST, tag=currenttag) # can also specify ids and idr to run just a subset of drout
