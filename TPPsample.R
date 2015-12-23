source("TPPmat.R")

library("deSolve")
library("stats") #contains optim()
library("lhs")
library("akima")
library("stringr")

taskids<-1:4
for (taskid in taskids)

currenttag <- paste0("20151223.",taskid)

Nsims_ds <- 500
Nsims_dr <- 200

targetepis <- list("Brazil"=c(52, 0.17, 0.014), "India"=c(195, 0.04, 0.022), "Philippines"=c(417, 0.002, 0.02), "SouthAfrica"=c(696, 0.61, 0.018)) # tb prev, HIV coprev, rrinc/inc

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)


values <- set.values()
mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
saveRDS(mergedvalues, file=paste0("genericvalues_",currenttag,".RDS"))

Nsamplepars_ds <- length(unlist(values$varied_ds))
Nsamplepars_dr <- length(unlist(values$varied_dr))

tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]

elementnames <- set.novelvalues()$elementnames

dsheader <- c("ids",  "targetprev","targetcoprev","targetdr",
              names(unlist(values)), 
            dssetup$statenames, tallynames) 
if(!file.exists(paste0("DScalibration_", currenttag, ".csv"))) { write(dsheader, sep =",", file=paste0("DScalibration_", currenttag, ".csv"), ncolumns=length(dsheader)) }

drheader <- c("ids","idr", "targetprev","targetcoprev","targetdr",
            names(unlist(values)), 
            drsetup$statenames, tallynames, paste0(tallynames,"10")) 
if(!file.exists(paste0("DRcalibration_", currenttag, ".csv"))) { write(drheader, sep = ",", file=paste0("DRcalibration_",currenttag,".csv"), ncolumns=length(drheader)) }

drtrajheader <- c("ids","idr", "targetprev","targetcoprev","targetdr",
              paste0( rep(-25:10, each=length(tallynames)), rep(tallynames, times=36) ))
if(!file.exists(paste0("DRtraj_", currenttag, ".csv"))) { write(drtrajheader, sep = ",", file=paste0("DRtraj_",currenttag,".csv"), ncolumns=length(drtrajheader)) }


LHS <- maximinLHS(Nsims_ds, Nsamplepars_ds); saveRDS(LHS, file=paste0("LHS_",currenttag,taskid,".RDS"))

for (isim in 1:Nsims_ds)
{
  dsvalues <- sample.values(values=values, whichparset="varied_ds", LHS=LHS, isim=isim)
  pars <- create.pars(setup = dssetup, values = dsvalues)
  
  # optimize for desired prev and coprev 
  optimat <- array(0, dim=c(0,4))
  b <- 4; h <- 0.0001; coprev <- 0; while(nrow(optimat)==0 | min(optimat[optimat[,2]==h/2, 4])< 0.5)
  {
    b <- max(4, ceiling(b/3)); prev <- 0; stopat4 <- F; while(prev<1000) 
    {
      dsvalues$cal$beta <- b; dsvalues$cal$hivrate <- h
      pars <- create.pars(dssetup, dsvalues)
      e <- equilib (pars=pars, tol=1)
      state <- e$log[nrow(e$log),2:(1+length(pars$fullpars$statenames))]
      prev <- sum(state[c(grep("^A", pars$fullpars$statenames), grep("^T", pars$fullpars$statenames))])
      coprev <- sum(state[c(grep("^A.+Hp", pars$fullpars$statenames), grep("^T.+Hp", pars$fullpars$statenames))]) / prev
      optimat <- rbind(optimat, c(b,h,prev,coprev))
      if (b==4 & prev>1000 & stopat4 ==F) { b <- 1; prev <- 0; stopat4 <- T} else if (b>14) b <- b+2 else if (prev<30) b <- b+2 else b <- b+1
    }
    print(paste0("Tried beta up to ", b, " for hivrate=", h))
    h <- h*2
    
  }                     
  
  #will treat each  epi (country) separately for the rest of the (DS and DR) calibration  
  for (tname in names(targetepis))
  {
    # find fit and create params for each set (country) of targetepis
    fit <- apply((t(optimat[,c(3,4)])-targetepis[[tname]][1:2])^2/targetepis[[tname]][1:2]^2, 2, sum)
    hs <- (0.0001*2^(0:10))[ 0.0001*2^(0:10) <= max(optimat[,2])]; bs <- min(optimat[,1]):max(optimat[,1])
  
    largefitmat <- interp(optimat[,1], 100*optimat[,2], z=fit, nx=50, ny=50, extrap=F, duplicate="mean")
    largefitmat$z[is.na(largefitmat$z)] <- 10000
    vindex <- which.min(largefitmat$z); aindex <- c(vindex - nrow(largefitmat$z)*floor(vindex/nrow(largefitmat$z)), ceiling(vindex/nrow(largefitmat$z)))
    
    dsvalues$cal$beta <- largefitmat$x[aindex[1]]; dsvalues$cal$hivrate <- largefitmat$y[aindex[2]]/100
    print(paste0("Chose beta=", dsvalues$cal$beta, ", hivrate=", dsvalues$cal$hivrate, " for sim #",isim))
    
    newpars <- create.pars(values=dsvalues, setup=dssetup)
    
    # and get equilibrium state
    opte <- equilib(pars=newpars, tol=0.1)
    dsstatenames <- dssetup$statenames
    estate <- with(opte,log[nrow(log),2:(length(dsstatenames)+1)])
  
    # save equilibrium state and values
    write(file=paste0("DScalibration_", currenttag, ".csv"), c(isim, unlist(targetepis[tname]), unlist(dsvalues), opte$log[nrow(opte$log),-1]), 
      sep=",", ncol=length(dsheader), append=TRUE)
    
  }
}
    
  
# add RifDR to newvalues, sampling again for acqres, transmissibility, and DSTrif
dsout <- read.csv(paste0("DScalibration_", currenttag, ".csv"), header=T)
for (isim in 1:Nsims_ds)
  for (tname in names(targetepis))
  {
    dsrow <- which(dsout$ids==isim & dsout$targetprev==(unlist(targetepis[tname])[1]) )
    estate <- dsout[dsrow , dssetup$statenames]
    drstate <- numeric(length(drsetup$statenames)); names(drstate) <- drsetup$statenames
    drstate[dsstatenames] <- estate
  
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
      drvalues <- sample.values(values=dsvalues, whichparset="varied_dr", LHS=drLHS, isim=isimdr)
      drpars <- create.pars(setup = drsetup, values = drvalues)
    
      # add and run to present over 25 years, with linear second line treatment scaleup to (present) max levels from years -10 to 0
      drend <- ode(unlist(drstate), seq(-25,10), dxdt, drpars$fullpars, rvary=T, nvary=F, do.tally=TRUE, method="adams")
        
      # will later screen out unacceptable range (or could adjust transmissibility to get targeted RifDR incidence)
  
      # save time-zero state and values
      results <- rbind(results, c(isim, isimdr, unlist(targetepis[tname]), 
                                  unlist(drvalues),
                                  drend[26,2:ncol(drend)],
                                  drend[36,tallynames]))
      write( c(isim, isimdr, unlist(targetepis[tname]), c(t(drend[, tallynames]))), file=paste0("DRtraj_",currenttag,".csv"), append=TRUE, sep=".", ncol=5+length(tallynames)*36 )
      print(paste0("Finished isimds=", isim, ", isimdr=", isimdr))
    } 

  write(t(results), ncolumns = ncol(results), append=TRUE,  sep = ",", file=paste0("DRcalibration_",currenttag,".csv"))
}


drout <- screendrout()


# limits output to those simulations with DR fraction of incidence in targetepi range

evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, targetpt="DS", DST="DSTall") # can also specify ids and idr to run just a subset of drout
evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, targetpt="DR", DST="DSTall", idr=1:10) # can also specify ids and idr to run just a subset of drout
evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, targetpt="DS", DST="DSTnone", idr=1:10) # can also specify ids and idr to run just a subset of drout
evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, targetpt="DR", DST="DSTnone", idr=1:10) # can also specify ids and idr to run just a subset of drout using ids and idr

