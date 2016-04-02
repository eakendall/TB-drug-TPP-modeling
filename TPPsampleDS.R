# In thi version, I sample an equally spaced grid for beta and log(hivrate), linearly interpolate between the sampled points to create a grid 5kx5k, then choose the best (sum of squared differences) fit on that 
# larger grid to the target prev and coprev.

taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1]
ntasks <- as.numeric(commandArgs(trailingOnly=TRUE))[2]
tname <- commandArgs(trailingOnly=TRUE)[3]
pessimistic <- commandArgs(trailingOnly=TRUE)[4]
location <- "../scratch/"

tag <- "20160313"
if (pessimistic) tag <- paste0(tag,"p")
Nsims_ds <- 500

source("TPPmat.R")
values <- set.values(pessimistic=pessimistic)

Nsamplepars_ds <- length(unlist(values$varied_ds)); ilimits <- ceiling(seq(0,Nsims_ds, length=ntasks+1)); print(ilimits)

if(file.exists(paste0(location,"LHS_",tag,".RDS"))) LHS <- readRDS(paste0(location,"LHS_",tag,".RDS")) else {LHS <- maximinLHS(Nsims_ds, Nsamplepars_ds); saveRDS(LHS, file=paste0(location,"LHS_",tag,".RDS"))}

if (Nsims_ds > nrow(LHS)) 
{ 
  oldLHS <- LHS; LHS <- augmentLHS(LHS, Nsims_ds - nrow(oldLHS)); saveRDS(LHS, file=paste0(location,"LHS_",tag,".RDS")); saveRDS(oldLHS, file=paste0(location,"oldLHS_",nrow(oldLHS),"_",tag,".RDS")) 
  ilimits <- ceiling(seq(nrow(oldLHS),Nsims_ds, length=ntasks+1)); print(ilimits)
}

currenttag <- paste0(tag,".",taskid)

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)

tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- set.novelvalues()$elementnames

dsheader <- c("ids",  "targetprev","targetcoprev","targetdr",  
              names(unlist(values)), 
            dssetup$statenames, tallynames) 
if(!file.exists(paste0(location,"DScalibration_", currenttag, ".csv"))) { write(dsheader, sep =",", file=paste0(location,"DScalibration_", currenttag, ".csv"), ncolumns=length(dsheader)) }

for (isim in (ilimits[taskid]+1):ilimits[taskid+1])
{
  dsvalues <- sample.values(values=values, whichparset="varied_ds", LHS=LHS, isim=isim, pessimistic=pessimistic)
  pars <- create.pars(setup = dssetup, values = dsvalues)
  
  # optimize for desired prev and coprev 
  optimat <- array(0, dim=c(0,4))
  h <- 0; coprev <- 0; while(h <= 0.0002 | min(optimat[optimat[,2]==h/4, 4])< 0.6)
  {
    b <- 2; prev <- 0; stopat3 <- F; while(prev<800) 
    {
      dsvalues$cal$beta <- b; dsvalues$cal$hivrate <- h
      pars <- create.pars(dssetup, dsvalues)
      e <- equilib(pars=pars, tol=1)
      state <- e$log[nrow(e$log),2:(1+length(pars$fullpars$statenames))]
      prev <- sum(state[c(grep("^A", pars$fullpars$statenames), grep("^T", pars$fullpars$statenames))])
      coprev <- sum(state[c(grep("^A.+Hp", pars$fullpars$statenames), grep("^T.+Hp", pars$fullpars$statenames))]) / prev
      optimat <- rbind(optimat, c(b,h,prev,coprev))
#       if (b==2 & prev>400 & stopat3 ==F) { b <- 1; prev <- 0; stopat3 <- T} else 
        if (b>9) b <- b+2 else 
          b <- b+1
    }
    print(paste0("Tried beta up to ", b, " for hivrate=", h))
    if (h==0) h <- 0.0001 else h <- h*2
  }                     
  
  saveRDS(optimat, file=paste0(location,"optimats",tag,"/optimat",isim,".RDS"))
  #will treat each  epi (country) separately for the rest of the (DS and DR) calibration  

optimat[,2][optimat[,2]==0] <- 0.000001
prevmat <- interp(optimat[,1], -log(optimat[,2]), z=optimat[,3], nx=5000, ny=5000, extrap=F, linear=TRUE, duplicate= "mean")
coprevmat <- interp(optimat[,1], -log(optimat[,2]), z=optimat[,4], nx=5000, ny=5000, extrap=F, linear=TRUE, duplicate= "mean")
fitmat <- (prevmat$z-targetepis[[tname]][1])^2/targetepis[[tname]][1]^2 + (coprevmat$z-targetepis[[tname]][2])^2/targetepis[[tname]][2]^2

vindex <- which.min(fitmat); aindex <- c(vindex - nrow(fitmat)*floor(vindex/nrow(fitmat)), ceiling(vindex/nrow(fitmat)))

dsvalues$cal$beta <- prevmat$x[aindex[1]]; dsvalues$cal$hivrate <- exp(-prevmat$y[aindex[2]])
if (dsvalues$cal$hivrate<0.00001) dsvalues$cal$hivrate <- 0

  print(paste0("Chose beta=", dsvalues$cal$beta, ", hivrate=", dsvalues$cal$hivrate, " for sim #",isim, " and epi of ", tname))
  
  newpars <- create.pars(values=dsvalues, setup=dssetup)
  
  # and get equilibrium state
  opte <- equilib(pars=newpars, tol=0.5)
  estate <- with(opte,log[nrow(log),2:(length(dssetup$statenames)+1)])
  
    # save equilibrium state and values
  write(file=paste0(location,"DScalibration_", currenttag, ".csv"), c(isim, unlist(targetepis[tname]), unlist(dsvalues), opte$log[nrow(opte$log),-1]), 
        sep=",", ncol=length(dsheader), append=TRUE)
  
#   }
  
}