taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1] #should have same array as dr runs
tname <- commandArgs(trailingOnly=TRUE)[2]
targetpt <- commandArgs(trailingOnly=TRUE)[3]
DST <- commandArgs(trailingOnly=TRUE)[4]

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

#idr <-  1:3 # use for test runs, or if splitting DR tasks into subtasks here

drout <- screendrout(drout_filename=paste0("fromMARCC/DRcalibration_",currenttag,".csv"), tolerance = 2)

# limits output to those simulations with DR fraction of incidence in targetepi range

evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, idr=idr, targetpt=targetpt, DST=DST, tag=currenttag) # can also specify ids and idr to run just a subset of drout
