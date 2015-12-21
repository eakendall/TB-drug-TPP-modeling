source("TPPmat.R")

library("deSolve")
library("stats") #contains optim()
library("lhs")
library("akima")
library("stringr")


currenttag <- "20151221"

Nsims_ds <- 3
Nsims_dr <- 3

targetepis <- list( "India"=c(195, 0.04, 0.022) )#list("Brazil"=c(52, 0.17, 0.014), "India"=c(195, 0.04, 0.022), "Philippines"=c(417, 0.002, 0.02), "SouthAfrica"=c(696, 0.61, 0.018)) # tb prev, HIV coprev, rrinc/inc

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)


values <- set.values()
mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
saveRDS(mergedvalues, file=paste0("genericvalues_",currenttag,".RDS"))

Nsamplepars_ds <- length(unlist(values$varied_ds))
Nsamplepars_dr <- length(unlist(values$varied_dr))

tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]

elementnames <- set.novelvalues$elementnames

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
    b <- max(4, b/3); prev <- 0; stopat4 <- F; while(prev<1000) 
    {
      dsvalues$cal$beta <- b; dsvalues$cal$hivrate <- h
      pars <- create.pars(dssetup, dsvalues)
      e <- equilib (pars=pars, tol=1)
      state <- e$log[nrow(e$log),2:(1+length(pars$fullpars$statenames))]
      prev <- sum(state[c(grep("^A", pars$fullpars$statenames), grep("^T", pars$fullpars$statenames))])
      coprev <- sum(state[c(grep("^A.+Hp", pars$fullpars$statenames), grep("^T.+Hp", pars$fullpars$statenames))]) / prev
      optimat <- rbind(optimat, c(b,h,prev,coprev))
      if (b==4 & prev>1000 & stopat4 ==F) { b <- 1; prev <- 0; stopat4 <- T} else if (b>20) b <- b+2 else if (prev<20) b <- b+2 else b <- b+1
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
#     fitmat <- array(0, dim=c(length(bs),length(hs))); dimnames(fitmat) <- list(bs, hs)
#     for (ifit in 1:nrow(optimat)) fitmat[as.character(optimat[ifit,1]), as.character(optimat[ifit,2])] <- fit[ifit]
#     fitmat[fitmat==0] <- max(fitmat)
#   
    # image(bs, hs, log(fitmat)) # useful for troubleshooting
  
    largefitmat <- interp(optimat[,1], 100*optimat[,2], z=fit, nx=50, ny=50, extrap=F, duplicate="mean")
    largefitmat$z[is.na(largefitmat$z)] <- 10000
    vindex <- which.min(largefitmat$z); aindex <- c(vindex - nrow(largefitmat$z)*floor(vindex/nrow(largefitmat$z)), ceiling(vindex/nrow(largefitmat$z)))
    
    dsvalues$cal$beta <- largefitmat$x[aindex[1]]; dsvalues$cal$hivrate <- largefitmat$y[aindex[2]]/100
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
      drend <- ode(drstate, seq(-25,0,by=0.1), dxdt, drpars$fullpars, rvary=T, nvary=F, do.tally=TRUE, method="adams")[251,]
        
      # will later screen out unacceptable range (or could adjust transmissibility to get targeted RifDR incidence)
  
      # save time-zero state and values
      results <- rbind(results, c(isim, isimdr, unlist(targetepis[tname]), 
                                  unlist(drvalues),
                                  drend[2:length(drend)]))
      print(paste0("Finished isimds=", isim, ", isimdr=", isimdr))
    } 

    write(t(results), ncolumns = ncol(results), append=TRUE,  sep = ",", file=paste0("DRcalibration_",currenttag,".csv"))
  }

}


drout <- screendrout()
extend_drout(drout)
dr10 <- cbind(drout, read.csv(paste0("DRextend_",currenttag,".csv"), header=T))
  
# limits output to those simulations with DR fraction of incidence in targetepi range

varynovel <- evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, targetpt="DS", DST="DSTall") # can also specify ids and idr to run just a subset of drout
# generates file of all TRP variations for all sampled simulations, for a given targetdt and dst, for each of potentially multiple targetepis

novelsubset <- loadnovel(targetpt="DS", DST="DSTall", targetepi=names(targetepis), tag=currenttag)
tail(novelsubset[,c(1:5, grep("tbdeaths", colnames(novelsubset)))]) #view

levels <- c("minimal","intermediate","optimal"); elementnames <- names(samplenovel(values))
traj <- array(0, dim=c(11,3,5)); dimnames(traj) <- list("t"=0:10, "level"=levels, "q"=c(0.025,0.25,0.5,0.075,0.975))
for (t in 0:10) for (l in levels) traj[t,l,] <- quantile(novelsubset[,novelheader==paste0(outcome, t, l)], c(0.025,0.25,0.5,0.075,0.975))

final_abs <- final_diff <- array(0,dim=c( length(elementnames)-1 , 2 , 5 )); dimnames(final_abs) <- dimnames(final_diff) <- list("vary"=names(genericTRP), "level"=c("minimal", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))
for (vary in elementnames[-length(elementnames)]) 
{
  final_abs[vary,1,] <- quantile(novelsubset[ , paste0(vary,".minimal.",outcome)], c(0.025,0.25,0.5,0.075,0.975))
  final_abs[vary,2,] <- quantile(novelsubset[ , paste0(vary,".optimal.",outcome)], c(0.025,0.25,0.5,0.075,0.975))
  

  final_diff[vary,1,] <- quantile(novelsubset[ , paste0("allintermediate10",outcome)] - novelsubset[ , paste0(vary,".",levelnames[1],".",outcome)]/
                                                                           novelsubset[ , paste0("allintermediate10",outcome)], c(0.025,0.25,0.5,0.075,0.975))
  final_diff[vary,2,] <- quantile(novelsubset[ , paste0("allintermediate10",outcome)] - novelsubset[ , paste0(vary,".",levelnames[2],".",outcome)]/
                                                                           novelsubset[ , paste0("allintermediate10",outcome)], c(0.025,0.25,0.5,0.075,0.975))
}  


}
        

plot(0:10, traj[,1,"0.5"])
points(0:10, traj[,2,"0.5"])
points(0:10, traj[,3,"0.5"])


barplot()

displayimpact <- function(drout, novelimpact)
{
  # plot: show basic median and 95% UR output before (e.g. from year -10) and after novel regimen x3 'all' levels: tb mort, tb incidence, panronsets, and time on rx x3 in same plot)
  # barplot: median tbmort for standard (none) TRP and for each minimal/optimal pair, 2x2 for DS/DR and w/wo DST
  # sensitivity analysis: PRCC of typical novel regimen impact (% tb mortality reduction) for each LHS-varied parameter;  
    ## and within targetpt DR / DST TRUE, PRCC for each TPP elements diff in tb mortality from minimal to optimal
  # resource use: barplots as above but for outcomes of diagnoses, DSTs (rif and novel in same plot), and rxmonths (all in same plot)

  
  inc <- drout[seq(1,nrow(drout), by=10,""]
  par(mfrow=c(2,2))  #plot baseline impact of novel regimen

