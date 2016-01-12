taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1]
ntasks <- as.numeric(commandArgs(trailingOnly=TRUE))[2]

tag <- "20160111"
Nsims_ds <- 250

source("TPPmat.R")
values <- set.values()

Nsamplepars_ds <- length(unlist(values$varied_ds)); #ilimits <- ceiling(seq(0,Nsims_ds, length=ntasks+1)); print(ilimits)

if(exists(paste0("LHS_",tag,".RDS"))) LHS <- readRDS(paste0("LHS_",tag,".RDS")) else 
  {LHS <- maximinLHS(Nsims_ds, Nsamplepars_ds); saveRDS(LHS, file=paste0("LHS_",tag,".RDS"))}

currenttag <- tag#paste0(tag,".",taskid)

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)

tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- set.novelvalues()$elementnames

dsheader <- c("ids",  "targetprev","targetcoprev","targetdr",  
              names(unlist(values)), 
            dssetup$statenames, tallynames) 
if(!file.exists(paste0("DScalibration_", currenttag, ".csv"))) { write(dsheader, sep =",", file=paste0("DScalibration_", currenttag, ".csv"), ncolumns=length(dsheader)) }

for (isim in 1:250)#(ilimits[taskid]+1):ilimits[taskid+1])
{
  dsvalues <- sample.values(values=values, whichparset="varied_ds", LHS=LHS, isim=isim)
  pars <- create.pars(setup = dssetup, values = dsvalues)
  
  # optimize for desired prev and coprev 
  optimat <- array(0, dim=c(0,4))
  h <- 0; coprev <- 0; while(h <= 0.0002 | min(optimat[optimat[,2]==h/4, 4])< 0.6)
  {
    b <- 3; prev <- 0; stopat3 <- F; while(prev<900) 
    {
      dsvalues$cal$beta <- b; dsvalues$cal$hivrate <- h
      pars <- create.pars(dssetup, dsvalues)
      e <- equilib(pars=pars, tol=1)
      state <- e$log[nrow(e$log),2:(1+length(pars$fullpars$statenames))]
      prev <- sum(state[c(grep("^A", pars$fullpars$statenames), grep("^T", pars$fullpars$statenames))])
      coprev <- sum(state[c(grep("^A.+Hp", pars$fullpars$statenames), grep("^T.+Hp", pars$fullpars$statenames))]) / prev
      optimat <- rbind(optimat, c(b,h,prev,coprev))
      if (b==3 & prev>400 & stopat3 ==F) { b <- 1; prev <- 0; stopat3 <- T} else 
        if (b>14) b <- b+2 else 
          b <- b+1
    }
    print(paste0("Tried beta up to ", b, " for hivrate=", h))
    if (h==0) h <- 0.0001 else h <- h*2
  }                     
  
  saveRDS(optimat, file=paste0("optimat",isim,".RDS"))
  #will treat each  epi (country) separately for the rest of the (DS and DR) calibration  

  for (tname in names(targetepis))
  {
    fit <- (optimat[,3]-targetepis[[tname]][1])^2/targetepis[[tname]][1]^2 + (optimat[,4]-targetepis[[tname]][2])^2/targetepis[[tname]][2]^2
    
    optimat[,2][optimat[,2]==0] <- 0.00000001
    
    largefitmat <- interp(optimat[,1], -log(optimat[,2]), z=fit, nx=50, ny=50, extrap=F, duplicate="mean")
    largefitmat$z[is.na(largefitmat$z)] <- 10000
    vindex <- which.min(largefitmat$z); aindex <- c(vindex - nrow(largefitmat$z)*floor(vindex/nrow(largefitmat$z)), ceiling(vindex/nrow(largefitmat$z)))
    
    dsvalues$cal$beta <- largefitmat$x[aindex[1]]; dsvalues$cal$hivrate <- exp(-largefitmat$y[aindex[2]])
    if (dsvalues$cal$hivrate<0.0000001) dsvalues$cal$hivrate <- 0
    
    print(paste0("Chose beta=", dsvalues$cal$beta, ", hivrate=", dsvalues$cal$hivrate, " for sim #",isim, " and epi of ", tname))
    
    newpars <- create.pars(values=dsvalues, setup=dssetup)
    
    # and get equilibrium state
    opte <- equilib(pars=newpars, tol=0.5)
    estate <- with(opte,log[nrow(log),2:(length(dssetup$statenames)+1)])
    
    # save equilibrium state and values
    write(file=paste0("DScalibration_", currenttag, ".csv"), c(isim, unlist(targetepis[tname]), unlist(dsvalues), opte$log[nrow(opte$log),-1]), 
          sep=",", ncol=length(dsheader), append=TRUE)
    
  }
  
}