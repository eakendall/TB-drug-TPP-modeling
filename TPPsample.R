library("stats") #contains optim()
library("lhs")
library("akima")




tag <- "20151212test2"

Nsims_ds <- 1
Nsims_dr <- 1

targetepis <- list( "India"=c(195, 0.04, 0.022) )#list("Brazil"=c(52, 0.17, 0.014), "India"=c(195, 0.04, 0.022), "Philippines"=c(417, 0.002, 0.02), "SouthAfrica"=c(696, 0.61, 0.018)) # tb prev, HIV coprev, rrinc/inc

Nsamplepars_ds <- nrow(samplepars("ds"))
Nsamplepars_dr <- nrow(samplepars("dr"))

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)

values <- set.values()
saveRDS(values, file=paste0("genericvalues_",tag,".RDS"))

tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]

dsheader <- c("ids",  "targetprev","targetcoprev","targetdr",
              "beta", "hivrate", names(unlist(values[unique(samplepars("ds")[,1])])), 
            dssetup$statenames) # tally names not yet defined, can pull from equilib colnames after time and dsstatenames
write(dsheader, sep =",", file=paste0("DScalibration_", tag, ".csv"), ncolumns=length(dsheader))

header <- c("ids","idr", "targetprev","targetcoprev","targetdr",
            "beta", "hivrate", names(unlist(values[unique(samplepars("ds")[,1])])), 
            names(unlist(values[unique(samplepars("dr")[,1])])), 
            drsetup$statenames, tallynames) # tally names not yet defined, can pull from equilib colnames after time and dsstatenames
write(header, sep = ",", file=paste0("DRcalibration_",tag,".csv"), ncolumns=length(header))



LHS <- maximinLHS(Nsims_ds, Nsamplepars_ds); saveRDS(LHS, file=paste0("LHS200_",tag,".RDS"))

for (isim in 1:Nsims_ds)
{
  dsvalues <- sample.values(values=values, whichparset="ds", LHS=LHS, isim=isim)
  pars <- create.pars(setup = dssetup, values = dsvalues)
  
  # optimize for desired prev and coprev 
  optimat <- array(0, dim=c(0,4))
  h <- 0.0001; coprev <- 0; while(nrow(optimat)==0 | min(optimat[optimat[,2]==h/2, 4])< 0.5)
  {
    b <- 4; prev <- 0; stopat4 <- F; while(prev<1000) 
    {
      dsvalues$beta <- b; dsvalues$hivrate <- h
      pars <- create.pars(dssetup, dsvalues)
      e <- equilib (pars=pars, tol=0.5)
      state <- e$log[nrow(e$log),2:(1+length(pars$fullpars$statenames))]
      prev <- sum(state[c(grep("^A", pars$fullpars$statenames), grep("^T", pars$fullpars$statenames))])
      coprev <- sum(state[c(grep("^A.+Hp", pars$fullpars$statenames), grep("^T.+Hp", pars$fullpars$statenames))]) / prev
      optimat <- rbind(optimat, c(b,h,prev,coprev))
      if (b==4 & prev>1000 & stopat4 ==F) { b <- 1; prev <- 0; stopat4 <- T} else if (b>20) b <- b+2 else if (prev<10) b <- b+2 else b <- b+1
    }
    h <- h*2
    print(paste0("Trying beta=", b, ", hivrate=", h))
  }                     
  
#will treat each  epi (country) separately for the rest of the (DS and DR) calibration  
  for (tname in names(targetepis))
  {
    # find fit and create params for each set (country) of targetepis
    
    fit <- apply((t(optimat[,c(3,4)])-targetepis[[tname]][1:2])^2/targetepis[[tname]][1:2]^2, 2, sum)
    hs <- (0.0001*2^(0:10))[ 0.0001*2^(0:10) <= max(optimat[,2])]; bs <- min(optimat[,1]):max(optimat[,1])
    fitmat <- array(0, dim=c(length(bs),length(hs))); dimnames(fitmat) <- list(bs, hs)
    for (ifit in 1:nrow(optimat)) fitmat[as.character(optimat[ifit,1]), as.character(optimat[ifit,2])] <- fit[ifit]
    fitmat[fitmat==0] <- max(fitmat)
  
    # image(bs, hs, log(fitmat)) # useful for troubleshooting
  
    largefitmat <- interp(optimat[,1], optimat[,2], z=fit, nx=50, ny=50, extrap=F, duplicate="mean")
    largefitmat$z[is.na(largefitmat$z)] <- 1000000
    vindex <- which.min(largefitmat$z); aindex <- c(vindex - nrow(largefitmat$z)*floor(vindex/nrow(largefitmat$z)), ceiling(vindex/nrow(largefitmat$z)))
    
    dsvalues$beta <- largefitmat$x[aindex[1]]; dsvalues$hivrate <- largefitmat$y[aindex[2]]
    print(paste0("Chose beta=", dsvalues$beta, ", hivrate=", dsvalues$hivrate))
    
    newpars <- create.pars(values=dsvalues, setup=dssetup)
    
    # and get equilibrium state
    opte <- equilib(func=dxdt, pars=newpars, tol=0.1)
    dsstatenames <- dssetup$statenames
    estate <- with(opte,log[nrow(log),2:(length(dsstatenames)+1)])
  
    # save equilibrium state and values
    write(file=paste0("DScalibration_", tag, ".csv"), c(isim, unlist(targetepis[tname]), dsvalues$beta, dsvalues$hivrate, unlist(dsvalues[unique(samplepars("ds")[,1])]), estate), 
      sep=",", ncol=4+length(unlist(values[unique(samplepars("ds")[,1])]))+2+length(dsstatenames), append=TRUE)
    
    
    # add RifDR to newvalues, sampling again for acqres, transmissibility, and DSTrif
  
    drstate <- numeric(length(drsetup$statenames)); names(drstate) <- drsetup$statenames
    drstate[dsstatenames] <- estate
  
    drLHS <- maximinLHS(Nsims_dr, Nsamplepars_dr)
    results <- numeric(0)
  
    for (isimdr in 1:Nsims_dr)
    {
      # sample 
      drvalues <- sample.values(values=dsvalues, whichparset="dr", LHS=drLHS, isim=isimdr)
      drpars <- create.pars(setup = drsetup, values = drvalues)
    
      # add and run to present over 25 years, with linear second line treatment scaleup to (present) max levels from years -5 to 0
      drend <- ode(drstate, seq(-25,0,by=0.1), dxdt, drpars$fullpars, rvary=T, nvary=F, do.tally=TRUE)[251,]
        
      # will later screen out unacceptable range (or could adjust transmissibility to get targeted RifDR incidence)
  
      # save time-zero state and values
      results <- rbind(results, c(isim, isimdr, unlist(targetepis[tname]), 
                                  newpars$values$beta, newpars$values$hivrate, unlist(newpars$values[unique(samplepars("ds")[,1])]), 
                                  unlist(newpars$values[unique(samplepars("dr")[,1])]), 
                                  drend[2:length(drend)]))
      print(paste0("Finished isimds=", isim, ", isimdr=", isimdr))
    }

    write(t(results), ncolumns = ncol(results), append=TRUE,  sep = ",", file=paste0("DRcalibration_",tag,".csv"))
  }

}



novelimpact <- evaltrp(genericvalues = values, drsetup = drsetup, drout=screendrout(), targetpt="DS", DST=FALSE)


  

