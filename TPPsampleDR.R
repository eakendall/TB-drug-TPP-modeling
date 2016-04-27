## !!!! Should select array to do only those (beyond 500) that aren't already done

taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1]
ntasks <- as.numeric(commandArgs(trailingOnly=TRUE))[2]
tname <- commandArgs(trailingOnly=TRUE)[3]
rDSTall <- commandArgs(trailingOnly=TRUE)[4]
Nsims_dr <- as.numeric(commandArgs(trailingOnly=TRUE))[5]
location <- "../scratch/"

tag <- "20160419p"
pessimistic <- TRUE

currenttag <- paste0(tname,"_",tag,".",taskid)
if (rDSTall==TRUE) currenttag <- paste0("rDSTall.",currenttag)

source("TPPmat.R")

if (tname != "India") dstag <- paste0(tname,"_",tag) else dstag <- tag
dsout <- list()
i <- 1; while(file.exists(paste0(location,"DScalibration_",dstag,".",i,".csv")))
{
  dsout <- rbind(dsout, read.csv(paste0(location,"DScalibration_",dstag,".",i,".csv"), header=TRUE))
  i <- i+1
}
dsout <- dsout[!duplicated(dsout$ids),] #want the earlier appearance of each ids
dsout <- dsout[order(dsout$ids),]

Nsims_ds <- max(dsout$ids); 
ilimits <- ceiling(seq(0,Nsims_ds, length=ntasks+1))


dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
values <- set.values(pessimistic=pessimistic)
mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
if (taskid==1) saveRDS(mergedvalues, file=paste0("genericvalues_",currenttag,".RDS"))

Nsamplepars_dr <- length(unlist(values$varied_dr))
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]

drheader <- c("ids","idr", "targetprev","targetcoprev","targetdr",
              names(unlist(values)), 
              drsetup$statenames, tallynames, paste0(tallynames,"10")) 
if(!file.exists(paste0(location,"DRcalibration_", currenttag, ".csv"))) { write(drheader, sep = ",", file=paste0(location,"DRcalibration_",currenttag,".csv"), ncolumns=length(drheader)) }

drtrajheader <- c("ids","idr", "targetprev","targetcoprev","targetdr",
                  paste0( rep(-25:10, each=length(tallynames)), rep(tallynames, times=36) ))
if(!file.exists(paste0(location,"DRtraj_", currenttag, ".csv"))) { write(drtrajheader, sep = ",", file=paste0(location,"DRtraj_",currenttag,".csv"), ncolumns=length(drtrajheader)) }

for (isim in (ilimits[taskid]+1):ilimits[taskid+1])
{
  dsrow <- which(dsout$ids==isim & dsout$targetprev==(unlist(targetepis[tname])[1]) )[1]
  estate <- dsout[dsrow , dssetup$statenames]
  drstate <- numeric(length(drsetup$statenames)); names(drstate) <- drsetup$statenames
  drstate[dssetup$statenames] <- estate
  
  dsvalues <- values
  valuevect <- dsout[dsrow, 4+(1:length(unlist(mergedvalues)))] 
  v <- 0; for (setname in names(dsvalues)) for (pname in names(dsvalues[[setname]]))
  { dsvalues[[setname]][[pname]] <- unlist(valuevect[v+(1:length(dsvalues[[setname]][[pname]]))]); #names(genericvalues[[pname]]) <- names()
    v <- v + length(dsvalues[[setname]][[pname]]) # move forward to start of next par vector in sampled values
  }    
  
  drLHS <- maximinLHS(Nsims_dr, Nsamplepars_dr)
  results <- numeric(0)
  
  for (isimdr in 1:Nsims_dr)
  {
    # sample 
    drvalues <- sample.values(values=dsvalues, whichparset="varied_dr", LHS=drLHS, isim=isimdr, pessimistic=pessimistic)
    if (rDSTall==TRUE) { drvalues$varied_dr$DSTrif_n <- 1; drvalues$varied_dr$noDSTrif_p <- 0 }
    drpars <- create.pars(setup = drsetup, values = drvalues)
    
    # add and run to present over 25 years, with linear second line treatment scaleup to (present) max levels from years -10 to 0
    drend <- ode(unlist(drstate), seq(-25,10), dxdt, drpars$fullpars, rvary=T, nvary=F, do.tally=TRUE, method="adams")
    
    # will later screen out unacceptable range (or could adjust transmissibility to get targeted RifDR incidence)
    
    # save time-zero state and values
    results <- rbind(results, c(isim, isimdr, unlist(targetepis[tname]), 
                                unlist(drvalues),
                                drend[26,2:ncol(drend)],
                                drend[36,tallynames]))
    write( c(isim, isimdr, unlist(targetepis[tname]), c(t(drend[, tallynames]))), file=paste0(location,"DRtraj_",currenttag,".csv"), append=TRUE, sep=",", ncol=5+length(tallynames)*36 )
    print(paste0("Finished isimds=", isim, ", isimdr=", isimdr))
  } 
  
  write(t(results), ncolumns = ncol(results), append=TRUE,  sep = ",", file=paste0(location,"DRcalibration_",currenttag,".csv"))
}
