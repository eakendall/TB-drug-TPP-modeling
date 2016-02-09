#source the function below, then run these:

taskid <- 1#as.numeric(commandArgs(trailingOnly=TRUE))[1]
tname <- "India"#commandArgs(trailingOnly=TRUE)[2]
targetpt <- "DS"#commandArgs(trailingOnly=TRUE)[3]
DST <- "DSTall"#commandArgs(trailingOnly=TRUE)[4]
if (targetpt=="DS") rDSTall <- TRUE else rDSTall <- FALSE

tag <- "20160201"
location<-""
dsout <- read.csv(paste0(location,"DScalibration_",tname,"_",tag,".1.csv"))
idrs <- 1:2

currenttag <- paste0(tname,"_",tag,".",taskid)
if(rDSTall) currenttag <- paste0("rDSTall.",currenttag)

source("TPPmat.R")
makemat <- function(pars) 
{
  with(pars, {
    mat <- array(0, dim=c(length(Tnames), length(Tnames), length(Rnames), length(Rnames), length(Hnames), length(Hnames))); dimnames(mat) <- list("1T"=Tnames, "2T"=Tnames, "1R"=Rnames, "2R"=Rnames, "1H"=Hnames, "2H"=Hnames)
    for (jr in Rnames) for (jt in Tnames) {  mat[jt, jt, jr, jr,"Hn","Hp"] <- mat[jt, jt, jr, jr,"Hn","Hp"] + hivrate }
    for (jr in Rnames) 
    {
      mat["An","S",jr,"R0","Hn","Hn"] <- mat["An","S",jr,"R0","Hn","Hn"] + selfcurerate;  # note all cures go back to R0 to simplify later computations
      for (jt in c("Ap","Ti")) { mat[jt,"C",jr,"R0","Hn","Hn"] <- mat[jt,"C",jr,"R0","Hn","Hn"] + selfcurerate } 
      for (jh in Hnames) {  
        mat["Ln","An",jr, jr, jh, jh] <- mat["Ln","An",jr, jr, jh, jh] + reactrate[jh]; 
        mat["Lp","Ap",jr, jr, jh, jh] <- mat["Lp","Ap",jr, jr, jh, jh] + reactrate[jh]  #reactivation, at hiv-dependent rate:
        diag(mat[,,jr,jr,jh,jh]) <- diag(mat[,,jr,jr,jh,jh]) + mort[jh]
        diag(mat[c("An","Ap","Ti"),c("An","Ap","Ti"),jr,jr,jh,jh]) <- diag(mat[c("An","Ap","Ti"),c("An","Ap","Ti"),jr,jr,jh,jh]) + tbmort[jh]
        mat["R","Ap",jr,jr,jh,jh] <- mat["R","Ap",jr,jr,jh,jh] + relapserate
    }  }
    
    relapsefracs <- function(period) # use relapse %s at 2, 4, and 6 months (defined in pars) for a 6 month regimen, and interpolate linearly 2--4 and 4--6 to get relapse % after any fraction of treatment course
    {
      fraction_completed <- (sum(Tperiod.lengths[(1:length(Tperiod.lengths))<=period]))/durations #changed to include whole period in which default occurs (will be 2nd period for all)
      relapse <- array(0, dim=c(length(Rnames),3)); dimnames(relapse) <- list(Rnames, regimens)
      relapse[,fraction_completed < 1/3] <- 1 
      relapse[,fraction_completed >= 1/3 & fraction_completed < 2/3] <- 
        (t((relapsemat) * (relapse246[2] + (relapse246[1]-relapse246[2])*(2/3-fraction_completed))))[, fraction_completed >= 1/3 & fraction_completed < 2/3]
      relapse[,fraction_completed >= 2/3 & fraction_completed < 1] <- 
        (t((relapsemat) * (relapse246[3] + (relapse246[2]-relapse246[3])*(1-fraction_completed))))[, fraction_completed >= 2/3 & fraction_completed < 1]
      relapse[,fraction_completed >= 1] <- t(relapsemat)[, fraction_completed >= 1]
      relapse[relapse>1] <- 1 #if sum exceeds 100%, set to 100% relapse
      return(relapse) 
    }
    
    ltfu <- ltfurate*durations*12
    for (jh in 1:length(Hnames))
    {
      if (length(usereg) > 0) for (jr in Rnames) for (reg in regimens[usereg])
      {
        periods <- paste0("T", reg, 1:Tperiods)
        jt <- 1
        while (sum(Tperiod.lengths[1:jt]) < durations[reg]-0.0001)
        {
          relapses <- relapsefracs(period=jt)
          mat[periods[jt],c(periods[jt+1],"R"),jr,jr,jh,jh] <- # move on; or discontinue and relapse 
            mat[periods[jt],c(periods[jt+1],"R"),jr,jr,jh,jh] + 1/Tperiod.lengths[jt]* c( (1-ltfu[reg]), (ltfu[reg])*relapses[jr, reg] )
          mat[periods[jt],"C",jr,"R0",jh,jh] <- # discontinue and cured (goes back to no resistance)
            mat[periods[jt],"C",jr,"R0",jh,jh] + 1/Tperiod.lengths[jt]* (ltfurate[reg])*(1-relapses[jr, reg])
          jt <- jt + 1
        }
        if (jt <= Tperiods) # run one more 
        {
          relapses <- relapsefracs(period=jt)
          mat[periods[jt],"R",jr,jr,jh,jh] <- 
            mat[periods[jt],"R",jr,jr,jh,jh] + 1/Tperiod.lengths[jt]*relapses[jr, reg]
          mat[periods[jt],"C",jr,"R0",jh,jh] <- 
            mat[periods[jt],"C",jr,"R0",jh,jh] + 1/Tperiod.lengths[jt]*(1-relapses[jr, reg])
        }
      }
    }
    return(mat)
  })
}


dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
values <- set.values()
mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))

tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]

drheader <- c("ids","idr", "targetprev","targetcoprev","targetdr",
              names(unlist(values)), 
              drsetup$statenames, tallynames, paste0(tallynames,"10")) 
if(!file.exists(paste0("DRcalibration_ltfu2mo.", currenttag, ".csv"))) { write(drheader, sep = ",", file=paste0("DRcalibration_ltfu2mo.",currenttag,".csv"), ncolumns=length(drheader)) }

drtrajheader <- c("ids","idr", "targetprev","targetcoprev","targetdr",
                  paste0( rep(-25:10, each=length(tallynames)), rep(tallynames, times=36) ))
if(!file.exists(paste0("DRtraj_ltfu2mo.", currenttag, ".csv"))) { write(drtrajheader, sep = ",", file=paste0("DRtraj_ltfu2mo.",currenttag,".csv"), ncolumns=length(drtrajheader)) }

olddrout <- read.csv(paste0("DRcalibration_",currenttag,".csv"))
  
isims <- min(olddrout$ids):max(olddrout$ids)

for (isim in isims)
{
  dsrow <- which(dsout$ids==isim & dsout$targetprev==(unlist(targetepis[tname])[1]) )[1]
  dsvalues <- values
  valuevect <- dsout[dsrow, 4+(1:length(unlist(mergedvalues)))] 
  v <- 0; for (setname in names(dsvalues)) for (pname in names(dsvalues[[setname]]))
  { dsvalues[[setname]][[pname]] <- unlist(valuevect[v+(1:length(dsvalues[[setname]][[pname]]))])
    v <- v + length(dsvalues[[setname]][[pname]]) # move forward to start of next par vector in sampled values
  }    
  pars <- create.pars(dssetup, dsvalues)
  e <- equilib(pars=pars, tol=1)
  estate <- e$log[nrow(e$log),2:(1+length(pars$fullpars$statenames))]
  
  drstate <- numeric(length(drsetup$statenames)); names(drstate) <- drsetup$statenames
  drstate[dssetup$statenames] <- estate
  
  Nsims_dr <-  max(olddrout$idr)
  
  for (isimdr in idrs)
  {
    drvalues <- values
    valuevect <- olddrout[which(olddrout$ids==isim & olddrout$idr==isimdr), 5+(1:length(unlist(mergedvalues)))] 
    v <- 0; for (setname in names(drvalues)) for (pname in names(drvalues[[setname]]))
    { drvalues[[setname]][[pname]] <- unlist(valuevect[v+(1:length(drvalues[[setname]][[pname]]))]); #names(genericvalues[[pname]]) <- names()
      v <- v + length(drvalues[[setname]][[pname]]) # move forward to start of next par vector in sampled values
    }    
    
    drpars <- create.pars(setup=drsetup, values=drvalues)
    
    drend <- ode(unlist(drstate), seq(-25,10), dxdt, drpars$fullpars, rvary=T, nvary=F, do.tally=TRUE, method="adams")
    
    write( c(isim, isimdr, unlist(targetepis[tname]), c(t(drend[, tallynames]))), file=paste0("DRtraj_ltfu2mo.",currenttag,".csv"), append=TRUE, sep=",", ncol=length(drtrajheader) )
    
    write(c(isim, isimdr, unlist(targetepis[tname]), 
            unlist(drvalues), drend[26,2:ncol(drend)],
            drend[36,tallynames]), append=TRUE,  sep = ",", ncol=length(drheader), file=paste0("DRcalibration_ltfu2mo.",currenttag,".csv"))
    
    print(paste0("Finished isimds=", isim, ", isimdr=", isimdr," for ltfu2mo ",currenttag))
  } 
  
}

drout <- screendrout(drout_filename = paste0("DRcalibration_ltfu2mo.",currenttag,".csv"), tolerance = 1.5)

evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, targetpt=targetpt, DST=DST, tag=paste0("ltfu2mo.",currenttag), location=location, rDSTall=rDSTall )

