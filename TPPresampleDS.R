tag <- "20160111"
Nsims_ds <- 250

source("TPPmat.R")

values <- set.values()

Nsamplepars_ds <- length(unlist(values$varied_ds))

LHS <- readRDS("LHS_20160105.RDS")
# if(exists(paste0("LHS_",tag,".RDS"))) LHS <- readRDS(paste0("LHS_",tag,".RDS")) else 
# {LHS <- maximinLHS(Nsims_ds, Nsamplepars_ds); saveRDS(LHS, file=paste0("LHS_",tag,".RDS"))}

currenttag <- tag

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)

tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- set.novelvalues()$elementnames

dsheader <- c("ids",  "targetprev","targetcoprev","targetdr",  
              names(unlist(values)), 
              dssetup$statenames, tallynames) 
if(!file.exists(paste0("DScalibration_", currenttag, ".csv"))) { write(dsheader, sep =",", file=paste0("DScalibration_", currenttag, ".csv"), ncolumns=length(dsheader)) }

for (isim in 1:Nsims_ds)
{
  dsvalues <- sample.values(values=values, whichparset="varied_ds", LHS=LHS, isim=isim)
  pars <- create.pars(setup = dssetup, values = dsvalues)
  
  optimat <- readRDS(paste0("optimats 20160106//optimat",isim,".RDS"))
  
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

