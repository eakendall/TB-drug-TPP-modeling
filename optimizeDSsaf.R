taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1]
ntasks <-as.numeric(commandArgs(trailingOnly=TRUE))[2]

tag <- "20160111"
Nsims_ds <- 250

source("TPPmat.R")
values <- set.values()

ilimits <- ceiling(seq(0,Nsims_ds, length=ntasks+1)); print(ilimits)

LHS <- readRDS(paste0("../scratch/LHS_20160105.RDS"))

currenttag <- paste0(tag,".",taskid)

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)

tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- set.novelvalues()$elementnames

targetarray <- array(unlist(targetepis), dim=c(3,4))

for (isim in (ilimits[taskid]+1):ilimits[taskid+1])
{
  dsvalues <- sample.values(values=values, whichparset="varied_ds", LHS=LHS, isim=isim)
  pars <- create.pars(setup = dssetup, values = dsvalues)
  
  # optimize for desired prev and coprev 
  optimat <- readRDS(paste0("../scratch/optimats/20160111/optimat",isim,".RDS"))
  h <- 0.256
  for (b in seq(2.5,5.5,by=1))
  {   dsvalues$cal$beta <- b; dsvalues$cal$hivrate <- h
      pars <- create.pars(dssetup, dsvalues)
      e <- equilib(pars=pars, tol=2)
      state <- e$log[nrow(e$log),2:(1+length(pars$fullpars$statenames))]
      prev <- sum(state[c(grep("^A", pars$fullpars$statenames), grep("^T", pars$fullpars$statenames))])
      coprev <- sum(state[c(grep("^A.+Hp", pars$fullpars$statenames), grep("^T.+Hp", pars$fullpars$statenames))]) / prev
      optimat <- rbind(optimat, c(b,h,prev,coprev))
   }
    print(paste0("Tried beta up to ", b, " for hivrate=", h))
  }                     
  
  saveRDS(optimat, file=paste0("../scratch/optimats/20160111/optimat",isim,".RDS"))
}