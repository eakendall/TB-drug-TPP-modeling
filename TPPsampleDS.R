taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1]
ntasks <- as.numeric(commandArgs(trailingOnly=TRUE))[2]

tag <- "20151223"
Nsims_ds <- 600


source("TPPmat.R")
values <- set.values()

Nsamplepars_ds <- length(unlist(values$varied_ds)); ilimits <- ceiling(seq(0,Nsims_ds, length=ntasks+1)); print(ilimits)

if(exists(paste0("LHS_",tag,".RDS"))) readRDS(paste0("LHS_",tag,".RDS")) else 
  {LHS <- maximinLHS(Nsims_ds, Nsamplepars_ds); saveRDS(LHS, file=paste0("LHS_",tag,".RDS"))}

currenttag <- paste0(tag,".",taskid)

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)

tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- set.novelvalues()$elementnames

dsheader <- c("ids",  "targetprev","targetcoprev","targetdr",
              names(unlist(values)), 
            dssetup$statenames, tallynames) 
if(!file.exists(paste0("DScalibration_", currenttag, ".csv"))) { write(dsheader, sep =",", file=paste0("../scratch/DScalibration_", currenttag, ".csv"), ncolumns=length(dsheader)) }

for (isim in (ilimits[taskid]+1):ilimits[taskid+1])
{
  dsvalues <- sample.values(values=values, whichparset="varied_ds", LHS=LHS, isim=isim)
  pars <- create.pars(setup = dssetup, values = dsvalues)
  
  # optimize for desired prev and coprev 
  optimat <- array(0, dim=c(0,4))
  b <- 4; h <- 0; coprev <- 0; while(min(optimat[optimat[,2]==h/2, 4])< 0.7)
  {
    b <- max(4, ceiling(b/3)); prev <- 0; stopat4 <- F; backtracked=F; while(prev<1000) 
    {
      dsvalues$cal$beta <- b; dsvalues$cal$hivrate <- h
      pars <- create.pars(dssetup, dsvalues)
      e <- equilib(pars=pars, tol=1)
      state <- e$log[nrow(e$log),2:(1+length(pars$fullpars$statenames))]
      prev <- sum(state[c(grep("^A", pars$fullpars$statenames), grep("^T", pars$fullpars$statenames))])
      coprev <- sum(state[c(grep("^A.+Hp", pars$fullpars$statenames), grep("^T.+Hp", pars$fullpars$statenames))]) / prev
      optimat <- rbind(optimat, c(b,h,prev,coprev))
      if (b==4 & prev>1000 & stopat4 ==F) { b <- 1; prev <- 0; stopat4 <- T} else 
        if (b==4 & prev>50 & coprev<0.2 & backtracked==F) {b <- 1; backtracked<-T} else
        if (b>14 | prev<30) b <- b+2 else 
          b <- b+1
    }
    print(paste0("Tried beta up to ", b, " for hivrate=", h))
    if (h==0) h <- 0.00025 else h <- h*2
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
    print(paste0("Chose beta=", dsvalues$cal$beta, ", hivrate=", dsvalues$cal$hivrate, " for sim #",isim, " and epi of ", tname))
    
    newpars <- create.pars(values=dsvalues, setup=dssetup)
    
    # and get equilibrium state
    opte <- equilib(pars=newpars, tol=0.1)
    dsstatenames <- dssetup$statenames
    estate <- with(opte,log[nrow(log),2:(length(dsstatenames)+1)])
  
    # save equilibrium state and values
    write(file=paste0("../scratch/DScalibration_", currenttag, ".csv"), c(isim, unlist(targetepis[tname]), unlist(dsvalues), opte$log[nrow(opte$log),-1]), 
      sep=",", ncol=length(dsheader), append=TRUE)
    
  }
}
    

