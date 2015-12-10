library("stats") #contains optim()
library("lhs")
library("akima")

# !! need to decide shape of prior distribution. and these numbers are just made up for now. 

samplepars <- function(whichparset)
{
  if (whichparset=="ds") samplepars <- t(array(c(
    "reactrate", 1, 0.0005,0.002,
    "rapidprog", 1, 0.05,0.3, 
    "tbmort",1, 0.1,0.3,
    "relapse_s",1,  0,0.1,
    "relapse246",2,  1.5, 5,
    "dxrate",1,  0.4,1.6    
  ), dim=c(4,6) )) else
  
  if (whichparset=="dr") samplepars <- t(array(c(
    "acqres_s", 1, 0.001, 0.009,
    "transmissibility", 1, 0.5, 0.9, 
    "DSTrif", 1, 0, 0.5,# how to handle treatment availability? !!
    "DSTrif", 2, 0.5, 1# how to handle treatment availability? !!
  ), dim=c(4, 4) )) else 
    
  stop("error, need to specify dr or ds")
  
  return(samplepars)
}


sample.values <- function(values, whichparset="ds", LHS, isim)
{
    if (missing(values)) values <- set.values()
    
  samplepars <-samplepars(whichparset)
  if (ncol(LHS) != nrow(samplepars)) stop("error, LHS and samplepar size mismatch")
    
  for (j in 1:ncol(LHS))
    values[[samplepars[j,1]]][as.numeric(samplepars[j,2])] <- as.numeric(samplepars[j,3]) + LHS[isim,j]*(as.numeric(samplepars[j,4])-as.numeric(samplepars[j,3]))
  
  return(values)
}



Nsims_ds <- 1
Nsims_dr <- 1

Nsamplepars_ds <- nrow(samplepars("ds"))
Nsamplepars_dr <- nrow(samplepars("dr"))

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)

values <- set.values()
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]

header <- c("ids","idr", 
            names(unlist(values[samplepars("ds")[,1]])), "beta", "hivrate", 
            names(unlist(values[samplepars("dr")[,1]])), 
            drsetup$statenames, tallynames) # tally names not yet defined, can pull from equilib colnames after time and dsstatenames
write(header, sep = ",", file="TPPcalibration.csv", ncolumns=length(header))

LHS <- maximinLHS(Nsims_ds, Nsamplepars_ds); save(LHS, file="LHS200_20151202.Rdata")

for (isim in 1:Nsims_ds)
{
  dsvalues <- sample.values(values=values, whichparset="ds", LHS=LHS, isim=isim)
  pars <- create.pars(setup = dssetup, values = dsvalues)
  
  # optimize for desired prev and coprev 
  optimat <- array(0, dim=c(0,4))
  h <- 0.0001; coprev <- 0; while(nrow(optimat)==0 | min(optimat[optimat[,2]==h/2, 4])<0.3)
  {
    b <- 4; prev <- 0; stopat4 <- F; while(prev<600) 
    {
      dsvalues$beta <- b; dsvalues$hivrate <- h
      pars <- create.pars(dssetup, dsvalues)
      e <- equilib (func=dxdt, pars=pars, tol=0.1)
      state <- e$log[nrow(e$log),2:(1+length(pars$fullpars$statenames))]
      prev <- sum(state[c(grep("^A", pars$fullpars$statenames), grep("^T", pars$fullpars$statenames))])
      coprev <- sum(state[c(grep("^A.+Hp", pars$fullpars$statenames), grep("^T.+Hp", pars$fullpars$statenames))]) / prev
      optimat <- rbind(optimat, c(b,h,prev,coprev))
      if (b==4 & prev>600 & stopat4 ==F) { b <- 1; prev <- 0; stopat4 <- T} else b <- b+1
    }
    h <- h*2
    print(paste0("Trying beta=", b, ", hivrate=", h))
  }                     
  
  targets <- c(200, 0.1) # tb prev, HIV coprev
  fit <- apply((t(optimat[,c(3,4)])-targets)^2/targets^2, 2, sum)
  hs <- (0.0001*2^(0:10))[ 0.0001*2^(0:10) <= max(optimat[,2])]; bs <- min(optimat[,1]):max(optimat[,1])
  fitmat <- array(0, dim=c(length(bs),length(hs))); dimnames(fitmat) <- list(bs, hs)
  for (ifit in 1:nrow(optimat)) fitmat[as.character(optimat[ifit,1]), as.character(optimat[ifit,2])] <- fit[ifit]
  fitmat[fitmat==0] <- max(fitmat)
  
#   image(bs, hs, log(fitmat)) # useful for troubleshooting
  
  largefitmat <- interp(optimat[,1], optimat[,2], z=fit, nx=50, ny=50, extrap=F, duplicate="mean")
  largefitmat$z[is.na(largefitmat$z)] <- 1000000
  vindex <- which.min(largefitmat$z); aindex <- c(vindex - nrow(largefitmat$z)*floor(vindex/nrow(largefitmat$z)), ceiling(vindex/nrow(largefitmat$z)))
  
  dsvalues$beta <- largefitmat$x[aindex[1]]; dsvalues$hivrate <- largefitmat$y[aindex[2]]
  print(paste0("Chose beta=", dsvalues$beta, ", hivrate=", dsvalues$hivrate))
    
  # and get equilibrium state
  newpars <- create.pars(values=dsvalues, setup=dssetup)
  opte <- equilib(func=dxdt, pars=newpars, tol=0.1)
  oldstatenames <- dssetup$statenames
  estate <- with(opte,log[nrow(log),2:(length(oldstatenames)+1)])
  
  # add MDR to newvalues, sampling again for acqres and transmissibility 
  
  newstate <- numeric(length(drsetup$statenames)); names(newstate) <- drsetup$statenames
  newstate[oldstatenames] <- estate
  
  drLHS <- maximinLHS(Nsims_dr, Nsamplepars_dr)
  results <- numeric(0)
  
  for (isimdr in 1:Nsims_dr)
  {
    # sample 
    newvalues <- sample.values(values=dsvalues, whichparset="dr", LHS=drLHS, isim=isimdr)
    newpars <- create.pars(setup = drsetup, values = newvalues)
    
    # add and run to present over 25 years
    mdrend <- advance(state=newstate, t0=0, dxdt, addedt=25, reportsteps = 1, fullpars=newpars$fullpars) 
    
    # !! will later screen out unacceptable range (or could adjust transmissibility to get targeted MDR incidence)
  
    # save time-zero state and values
    results <- rbind(results, c(isim, isimdr, 
                                unlist(newpars$values[samplepars("ds")[,1]]), newpars$values$beta, newpars$values$hivrate, 
                                unlist(newpars$values[samplepars("dr")[,1]]), 
                                mdrend[2,2:ncol(mdrend)]))
    print(paste0("Finished isimds=", isim, ", isimdr=", isimdr))
  }

  write(t(results), ncolumns = ncol(results), append=TRUE,  sep = ",", file="TPPcalibration.csv")
}
save(values, file="genericvalues_20151208")

screenmdrout <- function(mdrout_filename="TPPcalibration.csv", minratio, maxratio)
{
  mdrout <- read.csv(file = mdrout_filename, header = TRUE)
  screened <- mdrout[mdrout[,"rrinc"]/mdrout[,"inc"] > minratio & mdrout[,"rrinc"]/mdrout[,"inc"] < maxratio]  
  return(screened)
}

evaltrp(values, drsetup, screenmdrout(minratio=0.001, maxratio= 0.2))


evaltrp <- function(genericvalues, drsetup, mdrout_filename="TPPcalibration.csv", ids, idr, trpvalues)
{
  mdrout <- read.csv(file = mdrout_filename, header = TRUE)
  allsampledpars <-rbind(samplepars("ds"), samplepars("dr"))
  nsampledpars <- 0; for (name in allsampledpars[,1])  nsampledpars <- nsampledpars + length(genericvalues[[name]]) # count individual parameter elements
  
  novelsetup <- setup.model(DRera = TRUE, treatSL = TRUE, treatnovel = TRUE)
  
  if (missing(ids) | missing (idr)) rows <- 1:nrow(mdrout) else rows <- 1:nrow(mdrout)[mdrout[,"ids"] %in% ids & mdrout[,"idr"] %in% idr]
  for (i in rows)
  {
    iter <- mdrout[i,c("ids", "idr")] #will include this as part of returned output
    
    valuevect <- mdrout[i, 2+(1:nsampledpars)] 
    j <- 0 #initialize at start of valuevect
    for (name in allsampledpars[,1])
    {
      genericvalues[[name]] <- valuevect[i,j+(1:length(genericvalues[[name]]))]
      j <- j + length(genericvalues[[name]]) # move forward to start of next par vector in sampled values
    }    
    
    state <- mdrout[i, 2+nsampledpars + (1:length(drsetup$statenames))]
  
    # now everything is pulled back out from mdrout. Need to set up for including novel regimen. 
    newstate <- numeric(length(novelsetup$statenames)); names(newstate) <- novelsetup$statenames
    newstate[drsetup$statenames] <- as.list(state)
    
    
    # sample each TRP for the novel regimen
    !!
    
  
      
      
      ## implement novel regimen with a given TRP
      # assign companion resistance to specified fractions of 0 and r's
      s_cr <- nvalues$cres["rifs"]; r_dr <- nvalues$cres["rifr"]; novelstate <- newstate * rep( rep(c(1-s_cr, s_cr, 1-s_cr, s_cr, 1-r_cr, r_cr, 1-r_cr, r_cr), each=length(novelsetup$statenames)/16), 2)
      
      pars <- create.pars(setup = novelsetup, values = genericvalues)
      
    
    
    
    
  }
}
  
