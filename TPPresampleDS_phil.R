
# In this version, I sample an equally spaced grid for beta and log(hivrate), linearly interpolate between the sampled points to create a grid 5kx5k, then choose the best (sum of squared differences) fit on that 
# larger grid to the target prev and coprev.

taskid <- 1#as.numeric(commandArgs(trailingOnly=TRUE))[1]
ntasks <- 1#as.numeric(commandArgs(trailingOnly=TRUE))[2]
tname <- "Philippines"#commandArgs(trailingOnly=TRUE)[3]
location <- ""

tag <- "20160201"
currenttag <- paste0(tname,"_",tag,".",taskid)

Nsims_ds <- 250

source("TPPmat.R")

values <- set.values()

ilimits <- ceiling(seq(0,Nsims_ds, length=ntasks+1)); print(ilimits)

LHS <- readRDS("LHS_20160105.RDS")

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)

tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- set.novelvalues()$elementnames

dsheader <- c("ids",  "targetprev","targetcoprev","targetdr",  
              names(unlist(values)), 
              dssetup$statenames, tallynames) 
if(!file.exists(paste0(location,"DScalibration_", currenttag, ".csv"))) { write(dsheader, sep =",", file=paste0(location,"DScalibration_", currenttag, ".csv"), ncolumns=length(dsheader)) }

for (isim in (ilimits[taskid]+1):ilimits[taskid+1])
{
  dsvalues <- sample.values(values=values, whichparset="varied_ds", LHS=LHS, isim=isim)
  pars <- create.pars(setup = dssetup, values = dsvalues)
  
  optimat <- readRDS(paste0(location,"optimats 20160111/optimat",isim,".RDS")); colnames(optimat) <- c("beta","hivrate","prev","coprev")
  
  optimat <- optimat[optimat[,2]==0,]
  if (nrow(optimat)>5) optimat <- optimat[optimat[,3]>100,]
  
  optimat <- data.frame(optimat)
  fitbeta <- predict(lm(beta ~ prev, data=optimat), list(prev=targetepis[[tname]][1]))

  dsvalues$cal$beta <- fitbeta; dsvalues$cal$hivrate <- 0
  
  print(paste0("Chose beta=", dsvalues$cal$beta, ", hivrate=", dsvalues$cal$hivrate, " for sim #",isim, " and epi of ", tname))
  
  newpars <- create.pars(values=dsvalues, setup=dssetup)
  
  # and get equilibrium state
  opte <- equilib(pars=newpars, tol=1)
  estate <- with(opte,log[nrow(log),2:(length(dssetup$statenames)+1)])
  print(paste0("And got equilibrium prev of ",with(opte,log[nrow(log),"prev"])," and coprev of ",with(opte,log[nrow(log),"coprev"])))
  # save equilibrium state and values
  write(file=paste0(location,"DScalibration_", currenttag, ".csv"), c(isim, unlist(targetepis[tname]), unlist(dsvalues), opte$log[nrow(opte$log),-1]), 
        sep=",", ncol=length(dsheader), append=TRUE)
    
 }
