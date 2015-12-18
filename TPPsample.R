source("TPPmat.R")

library("stats") #contains optim()
library("lhs")
library("akima")
library("stringr")


currenttag <- "20151218"

Nsims_ds <- 3
Nsims_dr <- 3

targetepis <- list( "India"=c(195, 0.04, 0.022) )#list("Brazil"=c(52, 0.17, 0.014), "India"=c(195, 0.04, 0.022), "Philippines"=c(417, 0.002, 0.02), "SouthAfrica"=c(696, 0.61, 0.018)) # tb prev, HIV coprev, rrinc/inc

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)

values <- set.values()
mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
saveRDS(mergedvalues, file=paste0("genericvalues_",currenttag,".RDS"))

Nsamplepars_ds <- length(unlist(values$varied_ds))
Nsamplepars_dr <- length(unlist(values$varied_dr))

tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]

elementnames <- names(samplenovel(values)); levelnames <- names(samplenovel(values)[[1]])

dsheader <- c("ids",  "targetprev","targetcoprev","targetdr",
              names(unlist(values)), 
            dssetup$statenames, tallynames) 
if(!file.exists(paste0("DScalibration_", currenttag, ".csv"))) { write(dsheader, sep =",", file=paste0("DScalibration_", currenttag, ".csv"), ncolumns=length(dsheader)) }

drheader <- c("ids","idr", "targetprev","targetcoprev","targetdr",
            names(unlist(values)), 
            drsetup$statenames, tallynames) 
if(!file.exists(paste0("DRcalibration_", currenttag, ".csv"))) { write(drheader, sep = ",", file=paste0("DRcalibration_",currenttag,".csv"), ncolumns=length(drheader)) }



LHS <- maximinLHS(Nsims_ds, Nsamplepars_ds); saveRDS(LHS, file=paste0("LHS_",currenttag,".RDS"))

for (isim in 1:Nsims_ds)
{
  dsvalues <- sample.values(values=values, whichparset="varied_ds", LHS=LHS, isim=isim)
  pars <- create.pars(setup = dssetup, values = dsvalues)
  
  # optimize for desired prev and coprev 
  optimat <- array(0, dim=c(0,4))
  h <- 0.0001; coprev <- 0; while(nrow(optimat)==0 | min(optimat[optimat[,2]==h/2, 4])< 0.5)
  {
    b <- 4; prev <- 0; stopat4 <- F; while(prev<1000) 
    {
      dsvalues$cal$beta <- b; dsvalues$cal$hivrate <- h
      pars <- create.pars(dssetup, dsvalues)
      e <- equilib (pars=pars, tol=1)
      state <- e$log[nrow(e$log),2:(1+length(pars$fullpars$statenames))]
      prev <- sum(state[c(grep("^A", pars$fullpars$statenames), grep("^T", pars$fullpars$statenames))])
      coprev <- sum(state[c(grep("^A.+Hp", pars$fullpars$statenames), grep("^T.+Hp", pars$fullpars$statenames))]) / prev
      optimat <- rbind(optimat, c(b,h,prev,coprev))
      if (b==4 & prev>1000 & stopat4 ==F) { b <- 1; prev <- 0; stopat4 <- T} else if (b>20) b <- b+2 else if (prev<10) b <- b+2 else b <- b+1
    }
    print(paste0("Tried beta tp to ", b, " for hivrate=", h))
    h <- h*2
    
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
    
    dsvalues$cal$beta <- largefitmat$x[aindex[1]]; dsvalues$cal$hivrate <- largefitmat$y[aindex[2]]
    print(paste0("Chose beta=", dsvalues$cal$beta, ", hivrate=", dsvalues$cal$hivrate))
    
    newpars <- create.pars(values=dsvalues, setup=dssetup)
    
    # and get equilibrium state
    opte <- equilib(pars=newpars, tol=0.1)
    dsstatenames <- dssetup$statenames
    estate <- with(opte,log[nrow(log),2:(length(dsstatenames)+1)])
  
    # save equilibrium state and values
    write(file=paste0("DScalibration_", currenttag, ".csv"), c(isim, unlist(targetepis[tname]), unlist(dsvalues), opte$log[nrow(opte$log),-1]), 
      sep=",", ncol=length(dsheader), append=TRUE)
    
    
    # add RifDR to newvalues, sampling again for acqres, transmissibility, and DSTrif
  
    drstate <- numeric(length(drsetup$statenames)); names(drstate) <- drsetup$statenames
    drstate[dsstatenames] <- estate
  
    drLHS <- maximinLHS(Nsims_dr, Nsamplepars_dr)
    results <- numeric(0)
  
    for (isimdr in 1:Nsims_dr)
    {
      # sample 
      drvalues <- sample.values(values=dsvalues, whichparset="varied_dr", LHS=drLHS, isim=isimdr)
      drpars <- create.pars(setup = drsetup, values = drvalues)
    
      # add and run to present over 25 years, with linear second line treatment scaleup to (present) max levels from years -10 to 0
      drend <- ode(drstate, seq(-25,0,by=0.1), dxdt, drpars$fullpars, rvary=T, nvary=F, do.tally=TRUE, method=lsodes)[251,]
        
      # will later screen out unacceptable range (or could adjust transmissibility to get targeted RifDR incidence)
  
      # save time-zero state and values
      results <- rbind(results, c(isim, isimdr, unlist(targetepis[tname]), 
                                  unlist(newpars$values),
                                  drend[2:length(drend)]))
      print(paste0("Finished isimds=", isim, ", isimdr=", isimdr))
    } 

    write(t(results), ncolumns = ncol(results), append=TRUE,  sep = ",", file=paste0("DRcalibration_",currenttag,".csv"))
  }

}


drout <- screendrout()
# limits output to those simulations with DR fraction of incidence in targetepi range

varynovel <- evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, targetpt="DS", DST=TRUE) # can also specify ids and idr to run just a subset of drout
# generates file of all TRP variations for all sampled simulations, for a given targetdt and dst, for each of potentially multiple targetepis

loadnovel <- function(targetpt, DST, targetepi, tag=currenttag)
{
  novelout <- subset(read.csv(paste0("TRPoutput_", targetpt, DST,"_", tag,".csv")), targetprev==targetepis[[targetepi]][1] ) 
  return(novelout)
}

collateimpact <- function(targetpts=c("DS","DR"), DSTs=c(TRUE, FALSE), epis=names(targetepis), outcomes="tbdeaths", tag=currenttag)
{
  # for each targetpt, dst, varied.element, level, and targetepi, determine median, IQR, and 95% UR of all trajectories' TB mortality (can add other outcomes later)
  impact <- as.list(targetpts); names(impact) <- targetpts
  for (pt in targetpts)
  {
    impact[[pt]] <- as.list(DSTs); names(impact[[pt]]) <- DSTs
    for (DST in DSTs)
    {
      impact[[pt]][[DST]] <- as.list(names(targetepis)); names(impact[[pt]][[DST]]) <- names(targetepis)
      for (epi in names(targetepis))
      {
        novelsubset <- loadnovel(pt, DST, epi, tag)
        for (outcome in outcomes) #just tbmort for now
        {
          impact[[pt]][[DST]][[epi]][[outcome]] <- as.list(c("traj","final_abs","final_diff")); names(impact[[pt]][[DST]][[epi]][[outcome]]) <- (c("traj","final_abs","final_diff"))
          # trajectory with CIs for typical novel TRP
          impact[[pt]][[DST]][[epi]][[outcome]]$traj <- array(0, dim=c(5,11)); dimnames(traj) <- list("q"=c(0.025,0.25,0.5,0.075,0.975), "t"=0:10)
          for (t in 0:10) impact[[pt]][[DST]][[epi]][[outcome]]$traj[,t] <- quantile(novelsubset[,paste0(outcome,t)], c(0.025,0.25,0.5,0.075,0.975))

          # final result with CIs for each TRP variation: as an absolute outcome, and as compared to the none TRP for the same values
          impact[[pt]][[DST]][[epi]][[outcome]]$final_abs <- impact[[pt]][[DST]][[epi]][[outcome]]$final_diff <- array(0,dim=c(5,length(elementnames),2))
          dimnames(impact[[pt]][[DST]][[epi]][[outcome]]$final_abs) <- dimnames(impact[[pt]][[DST]][[epi]][[outcome]]$final_diff) <-                                                                             final_diff) <- 
            list("q"=c(0.025,0.25,0.5,0.075,0.975), "vary"=elementnames, "level"=levelnames)
          for (vary in elementnames) 
          {
            impact[[pt]][[DST]][[epi]][[outcome]]$final_abs[,vary,1] <- quantile(novelsubset[ ,paste0(vary,".",levelnames[1],".",outcome)], c(0.025,0.25,0.5,0.075,0.975))
            impact[[pt]][[DST]][[epi]][[outcome]]$final_abs[,vary,2] <- quantile(novelsubset[ ,paste0(vary,".",levelnames[2],".",outcome)], c(0.025,0.25,0.5,0.075,0.975))

            impact[[pt]][[DST]][[epi]][[outcome]]$final_diff[,vary,1] <- quantile( novelsubset[ ,paste0("none.none.",outcome)] - novelsubset[ ,paste0(vary,".",levelnames[1],".",outcome)]/
                                               novelsubset[ ,paste0("none.none.",outcome)], c(0.025,0.25,0.5,0.075,0.975))
            impact[[pt]][[DST]][[epi]][[outcome]]$final_diff[,vary,2] <- quantile( novelsubset[ ,paste0("none.none.",outcome)] - novelsubset[ ,paste0(vary,".",levelnames[2],".",outcome)]/
                                               novelsubset[ ,paste0("none.none.",outcome)], c(0.025,0.25,0.5,0.075,0.975))
          }
        }
      }
    }
  }
  return(impact)
}
        

displayimpact <- function(drout, novelimpact)
{
  # plot: show basic median and 95% UR output before (e.g. from year -10) and after novel regimen introduction: tb mort, tb incidence, panronsets, and time on rx x3 in same plot)
  # barplot: median tbmort for standard (none) TRP and for each minimal/optimal pair, 2x2 for DS/DR and w/wo DST
  # sensitivity analysis: PRCC of typical novel regimen impact (% tb mortality reduction) for each LHS-varied parameter;  
    ## and within targetpt DR / DST TRUE, PRCC for each TPP elements diff in tb mortality from minimal to optimal
  # resource use: barplots as above but for outcomes of diagnoses, DSTs (rif and novel in same plot), and rxmonths (all in same plot)

  
  inc <- drout[seq(1,nrow(drout), by=10,""]
  par(mfrow=c(2,2))  #plot baseline impact of novel regimen

