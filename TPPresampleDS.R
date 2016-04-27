taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1]
ntasks <- as.numeric(commandArgs(trailingOnly=TRUE))[2]
tname <- commandArgs(trailingOnly=TRUE)[3]
location <- "../scratch/"

tag <- "20160419p"
pessimistic <- TRUE
currenttag <- paste0(tname,"_",tag,".",taskid)

# LHS <- readRDS(paste0("LHS_",tag,".RDS"))
# LHS <- readRDS(paste0("oldLHS_100_",tag,".RDS"))
Nsims_ds <- 100

source("TPPmat.R")

values <- set.values()
genericvalues <- values

ilimits <- ceiling(seq(0,Nsims_ds, length=ntasks+1)); print(ilimits)

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)

tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- set.novelvalues()$elementnames

dsheader <- c("ids",  "targetprev","targetcoprev","targetdr",  
              names(unlist(values)), 
              dssetup$statenames, tallynames) 
if(!file.exists(paste0(location,"DScalibration_", currenttag, ".csv"))) { write(dsheader, sep =",", file=paste0(location,"DScalibration_", currenttag, ".csv"), ncolumns=length(dsheader)) }

dsout <- list(); i <- 1; while(file.exists(paste0("DScalibration_20160419p.",i,".csv"))) {dsout <- rbind(dsout, read.csv(paste0("DScalibration_20160419p.",i,".csv"))); i <- i+1}

for (isim in (ilimits[taskid]+1):ilimits[taskid+1])
{
  valuevect <- dsout[which(dsout$ids==isim),names(unlist(genericvalues))]
  dsvalues <- genericvalues
  v <- 0; for (setname in names(genericvalues)) for (pname in names(genericvalues[[setname]]))
  { 
    dsvalues[[setname]][[pname]] <- unlist(valuevect[v+(1:length(dsvalues[[setname]][[pname]]))]); #names(genericvalues[[pname]]) <- names()
    v <- v + length(dsvalues[[setname]][[pname]]) # move forward to start of next par vector in sampled values
  }    
  
#   dsvalues <- 
#     sample.values(values=values, whichparset="varied_ds", LHS=LHS, isim=isim, pessimistic=pessimistic)
  pars <- create.pars(setup = dssetup, values = dsvalues)
  
  optimat <- readRDS(paste0(location,"optimats",tag,"/optimat",isim,".RDS"))
  
  if (tname == "Philippines")
  {  
    usethese <- optimat[optimat[,2]==0&optimat[,3]>100,]
    colnames(usethese) <- c("b","h","prev","coprev")
    model <- lm(b ~ prev, data=data.frame(usethese))
    dsvalues$cal$hivrate <- 0; dsvalues$cal$beta <- as.numeric(model$coefficients[1] + model$coefficients[2]*targetepis[[tname]][1])
  } else
  if (tname != "Philippines")
  { optimat[,2][optimat[,2]==0] <- 0.000001
    prevmat <- interp(optimat[,1], -log(optimat[,2]), z=optimat[,3], nx=5000, ny=5000, extrap=F, linear=TRUE, duplicate= "mean")
    coprevmat <- interp(optimat[,1], -log(optimat[,2]), z=optimat[,4], nx=5000, ny=5000, extrap=F, linear=TRUE, duplicate= "mean")
    fitmat <- (prevmat$z-targetepis[[tname]][1])^2/targetepis[[tname]][1]^2 + (coprevmat$z-targetepis[[tname]][2])^2/targetepis[[tname]][2]^2
  
    vindex <- which.min(fitmat); aindex <- c(vindex - nrow(fitmat)*floor(vindex/nrow(fitmat)), ceiling(vindex/nrow(fitmat)))
  
    dsvalues$cal$beta <- prevmat$x[aindex[1]]; dsvalues$cal$hivrate <- exp(-prevmat$y[aindex[2]])
  }
  
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
