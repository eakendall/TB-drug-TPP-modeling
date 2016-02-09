# Want to track resistance-related outcomes over time: aprev=excluding those suppressed on treatment, rprev, cprev, nprev, novelprev, rnovelprev, onsets, ronsets, nronsets, tbdeaths, rrdeaths
# for an otherwise-intermediate regimen with all combinations (3x3) of companion and barrier
# without and with novel DST, for India
# and maybe for DS and DR regimen (but will do just DS first)


taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1] #the idr's we want to run
tname <- commandArgs(trailingOnly=TRUE)[2]
targetpt <- commandArgs(trailingOnly=TRUE)[3]
DST <- commandArgs(trailingOnly=TRUE)[4]
rDSTall <- commandArgs(trailingOnly=TRUE)[5]
location<-"../scratch/"

tag <- "20160201"
currenttag <- paste0(tname,"_",tag)
if (rDSTall == TRUE) currenttag <- paste0("rDSTall.",currenttag)
tasktag <- paste0(currenttag,".idr",taskid)
source("TPPmat.R")

dxdt <- function(t, state, fullpars, rvary, nvary, do.tally=FALSE)
{
  DSTrif_t <- numeric(2); availability_t <- numeric(1)
  
  if (missing(rvary)) { if (2 %in% fullpars$usereg) {rvary<-TRUE} else {rvary<-FALSE} }
  if (rvary==TRUE) { if (t <= -10) {DSTrif_t <- 0*fullpars$DSTrif} else if (t <= 0 && t > -10) {DSTrif_t <- (10+t)/10 * fullpars$DSTrif} else if (t > 0) {DSTrif_t <- fullpars$DSTrif}
  } else {DSTrif_t <- fullpars$DSTrif}
  
  if (missing(nvary)) { if (3 %in% fullpars$usereg) {nvary<-TRUE} else {nvary<-FALSE} }
  if (nvary==TRUE) { if (t <= 0) availability_t <- 0 else if (t <= 3 && t > 0) availability_t <- (t/3)*fullpars$availability else if (t > 3) availability_t <- fullpars$availability 
  } else availability_t <- fullpars$availability # 3 here is the scale-up time
  
  with(fullpars, {
    if (length(mat) != length(state)^2) {stop("Error: Initial-state and transition-matrix size mismatch.")}
    
    statemat <- array(unlist(state), dim=c(length(Tnames), length(Rnames), length(Hnames))); dimnames(statemat) <- list(Tnames, Rnames, Hnames)
    if (length(Rnames)==1) { FOI <- sum(statemat[c("An","Ap","Ti"),,]) * transmissibility * beta / sum(statemat) 
    } else  { FOI <- apply(statemat[c("An","Ap","Ti"),,], 2, sum) * transmissibility * beta / sum(statemat) }# FOI by strain
    names(FOI) <- names(transmissibility) <- Rnames
    
    # infection
    #      if (length(Rnames)>1 && sum(statemat[c("S","C"), -1, ]) > 0.0001 ) { stop("Error: Some susceptibles have drug resistance and won't be included in infection events.") }
    for (jh in Hnames)
    {
      mat["S","Ln","R0",,jh,jh] <- mat["S","Ln","R0",,jh,jh] + FOI*(1-rapidprog[jh])
      mat["S","An","R0",,jh,jh] <- mat["S","An","R0",,jh,jh] + FOI*rapidprog[jh]
      mat["C","Lp","R0",,jh,jh] <- mat["C","Lp","R0",,jh,jh] + FOI*(1-rapidprog[jh])
      mat["C","Ap","R0",,jh,jh] <- mat["C","Ap","R0",,jh,jh] + FOI*rapidprog[jh]
      
      # superinfection
      for (jr in Rnames) # unlke above, vector now refers to resistance phenotype of *origin*
      {
        mat["Ln","Ln",,jr,jh,jh] <- mat["Ln","Ln",,jr,jh,jh] + FOI[jr]*(1-rapidprog[jh]*latreduct)*transmissibility[jr]/(transmissibility + transmissibility[jr])
        mat["Ln","An",,jr,jh,jh] <- mat["Ln","An",,jr,jh,jh] + FOI[jr]*rapidprog[jh]*latreduct
      }
    }
    
    ## TB diagnosis and treatment initiation
    for (jh in 1:length(Hnames))
    {
      # tb is diagnosed and rif DST is performed depending on An vs Ap, 
      # then regimen is chosen according to RifDR diagnosis, new regimen's target population, new regimen DST use, and medical eligibility by HIV status
      # considers all resistance including those not included in model setup
      startmat <- array(0,dim=c(2,8,3)); dimnames(startmat) <- list("thist"=c("An","Ap"), c("R0", "Rc", "Rn", "Rcn", "Rr", "Rrc", "Rrn", "Rrcn"), regimens) # probability of starting new regimen, by n/p treamtent history (down) and resistance (across) 
      
      startmat[,,"n"] <- cbind( 
        outer( targetpop[1] * c(1,1) * availability_t * eligibility[jh] , c(1, 1-DSTnew[1], 1-DSTnew[2], 1-max(DSTnew)) ) , # fraction of non-RifDR (Rnames 1:4) that get novel regimen
        outer( (targetpop[1]*(1-DSTrif_t) + targetpop[2]*DSTrif_t) * availability_t * eligibility[jh] , c(1, 1-DSTnew[1], 1-DSTnew[2], 1-max(DSTnew)) ) # fraction of undiagnosed RifDR can get new DS regimen, or diagnosed RifDR can get new RifDR regimen
      ) 
      startmat[,1:4,"s"] <- 1 - startmat[,1:4,"n"]; startmat[,5:8,"s"] <- (1-DSTrif_t) * (1 - startmat[,5:8,"n"]) 
      startmat[,5:8,"r"] <- DSTrif_t * (1 - startmat[,5:8,"n"])
      
      # now assign to ineffective treatment, or pending relapse with acquired resistance (of those who don't fail outright), or to first phase of effective treatment
      for (it in c("An", "Ap")) 
      {
        for (nreg in usereg)
        {
          startrate <- (1-initialloss[nreg])*dxrate[it, jh] #initial loss to follow up remains active; this allows it to differ by regimen (or allows us to turn off a regimen e.g. have no alternative to the novel regimen)
          #acquired resistance moves to pending relapse for *new* strain
          mat[it, "R", , , jh, jh]  <- 
            mat[it, "R", , , jh, jh] + startrate * startmat[it,Rnames,nreg] * acqresmat[,,nreg]
          
          #of those who don't acquire resistance (1-rowSums(acqresmat)), will split between Ti (failmat) and T1 for the selected (startmat) regimen
          if (length(Rnames)>1)
          {
            diag(mat[it, "Ti", , , jh, jh]) <- 
              diag(mat[it, "Ti", , , jh, jh]) + startrate * startmat[it, Rnames ,nreg] * (1-rowSums(acqresmat[,,nreg]))*failmat[nreg, ] #failures go to ineffective treatment (with ongoing infectiousness and increased mortality risk)
            diag(mat[it, grep("T[srn]1",Tnames)[nreg], , , jh, jh]) <- 
              diag(mat[it, grep("T[srn]1",Tnames)[nreg], , , jh, jh]) + startrate * startmat[it, Rnames, nreg] * (1-rowSums(acqresmat[,,nreg]))*(1 - failmat[nreg, ]) #for the remainder, treatment is initially effective (i.e. outcome will be cure vs relapse if they complete at least 2 months)
          } else #if only the standard regimen is an option, diag and rowSums functions cause errors
          {
            mat[it, "Ti", 1,1, jh, jh] <- 
              mat[it, "Ti", 1,1, jh, jh] + startrate * startmat[it, Rnames ,nreg] * (1-(acqresmat[,,nreg]))*failmat[nreg,Rnames ] #failures go to ineffective treatment (with ongoing infectiousness and increased mortality risk)
            mat[it, grep("T[srn]1",Tnames)[nreg], , , jh, jh] <- 
              mat[it, grep("T[srn]1",Tnames)[nreg], , , jh, jh] +  startrate * startmat[it, Rnames, nreg] * (1-(acqresmat[,,nreg]))*(1 - failmat[nreg, Rnames]) #for the remainder, treatment is initially effective (i.e. outcome will be cure vs relapse if they complete at least 2 months)            
          }
        }
      }
    }
    
    
    ########### transform mat
    newmat <- aperm(mat, c(1,3,5,2,4,6)) # this will translate into a 2d matrix more easily
    squaremat <- array(newmat, dim=c(length(Tnames)*length(Rnames)*length(Hnames), length(Tnames)*length(Rnames)*length(Hnames))) #2d state1 to state2 transmition matrix
    dimnames(squaremat) <- list(statenames, statenames)
    
    
    ########## tally how much each current state contributes to outcomes of interest (then will multiply by state and append to output)
    
    outcomes <- c("aprev", "rprev", "cprev", "nprev", "novelprev", "rnovelprev", "onsets", "ronsets", "nonsets", "tbdeaths", "rrdeaths")
    tally <- array(0,dim=c(length(Tnames)*length(Rnames)*length(Hnames), length(outcomes))); dimnames(tally) <- list(statenames, outcomes)
    
    if (do.tally==TRUE)
    {
      tally[c(grep("^A", statenames), grep("^Ti", statenames)), "aprev"] <- 1
      
      tally[c(grep("^A.[.]Rr", statenames), grep("^Ti[.]Rr", statenames)), "rprev"] <- 1
      
      tally[c(grep("^A.[.]R(rc|c)", statenames), grep("^Ti[.]R(rc|c)", statenames)), "cprev"] <- 1
      
      tally[c(grep("^A.[.]R(rcn|cn|rn|n)", statenames), grep("^Ti[.]R(rcn|cn|rn|n)", statenames)), "nprev"] <- 1
      
      tally[c(grep("^A.[.]R(rc|rn|c|n)", statenames), grep("^Ti[.]R(rc|rn|c|n)", statenames)), "novelprev"] <- 1
      
      tally[c(grep("^A.[.]R(rc|rn)", statenames), grep("^Ti[.]R(rc|rn)", statenames)), "rnovelprev"] <- 1
      
      tally[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames),grep("^R",statenames)),"onsets"] <- #doens't include relapses (mostly early), but does include reinfections (mostly late) - i.3. this is a count of new infecitons, whether or not they are recognized as such
        apply(squaremat[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames),grep("^R",statenames)),  grep("^A",statenames)], 1, sum)
      
      tally[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames), grep("^R",statenames)),"ronsets"] <- #this alternative also includes resistance relapses (including resistance acquisitions with relapse)
        apply(squaremat[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames),grep("^R",statenames)),  grep("^A.[.]Rr",statenames)], 1, sum)
      
      tally[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames), grep("^R",statenames)),"nonsets"] <-
        apply(squaremat[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames),grep("^R",statenames)),  grep("^A.[.]R(rcn|cn|rn|n)",statenames)], 1, sum)
      
      tally[c(grep("^A", statenames), grep("^Ti", statenames)), "tbdeaths"] <- rep(tbmort, each=(length(c(grep("^A", statenames), grep("^Ti", statenames)))/2))
      
      tally[c(grep("^A.[.]Rr", statenames), grep("^Ti[.]Rr", statenames)), "rrdeaths"] <- rep(tbmort, each=(length(c(grep("^A.[.]Rr", statenames), grep("^Ti[.]Rr", statenames)))/2))
    }
    
    tallied <- t(tally) %*% state ; names(tallied) <- outcomes; 
    
    ## combine state vector and change matrix to get dxdt
    dxdt <- numeric(length(state))
    deaths <- state*c(diag(squaremat)); diag(squaremat) <- 0 #store deaths elsewhere
    dxdt <- dxdt + t(squaremat) %*% state - rowSums(squaremat*state) - deaths #ins minus outs minus deaths
    
    #   births:
    # assume relatively flat total and dr tb foi over the past 15 years (as this will be true by time 0 apart from <50% overestimating DR exposure), and 
    # estimate as an exponential the latent TB in 15yos enterting the population. Assume no latent novel drug resistance. 
    # And companion drug resistance will be appropriately assigned to each fraction at time 0, but after that it will need to continue: 
    # so keep the same ratios of c and not c among new latents after time 0, as were assigned at time zero.
    
    
    if (length(Rnames)==1)    
    {
      dxdt[statenames=="S.R0.Hn"] <- dxdt[statenames=="S.R0.Hn"] + sum(deaths) * exp(-15*sum(FOI))
      dxdt[statenames=="Ln.R0.Hn"] <- dxdt[statenames=="Ln.R0.Hn"] + sum(deaths) * (1-exp(-15*sum(FOI)))
    }
    
    if (length(Rnames)==2 | (t<=0  & length(Rnames)==8))
    {
      dxdt[statenames=="S.R0.Hn"] <- dxdt[statenames=="S.R0.Hn"] + sum(deaths) * exp(-15*sum(FOI))
      dxdt[statenames %in% c("Ln.R0.Hn", "Ln.Rr.Hn")] <- dxdt[statenames %in% c("Ln.R0.Hn", "Ln.Rr.Hn")] + sum(deaths) * 
        c( sum(FOI[grep("R[0cn]+", names(FOI))]), sum(FOI[grep("Rr+", names(FOI))]) )/ sum(FOI) * (1-exp(-15*sum(FOI)))
    }
    
    
    if (t>0 & length(Rnames)==8) 
    {
      dxdt[statenames=="S.R0.Hn"] <- dxdt[statenames=="S.R0.Hn"] + sum(deaths) * exp(-15*sum(FOI))
      dxdt[statenames %in% c("Ln.R0.Hn", "Ln.Rr.Hn", "Ln.Rc.Hn", "Ln.Rrc.Hn")] <- dxdt[statenames %in% c("Ln.R0.Hn", "Ln.Rr.Hn", "Ln.Rc.Hn", "Ln.Rrc.Hn")] + sum(deaths) * 
        c( sum(FOI[c(1,3)]), sum(FOI[c(5,7)]), sum(FOI[c(2,4)]), sum(FOI[c(6,8)]) )/ sum(FOI) * c(1-fullpars$cres, fullpars$cres) * (1-exp(-15*sum(FOI)))
    }
    
    # return dxdt to ode, and also return tally for tracking purposes
    return(list(dxdt, tallied))
  }
  )
}


dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
values <- set.values()
genericvalues <- mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
rtallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- c("all",set.novelvalues()$elementnames)

alldrout <- numeric(0)
i <- 1; while(file.exists(paste0(location,"DRcalibration_",currenttag,".",i,".csv")))
{alldrout <- rbind(alldrout, read.csv(paste0(location,"DRcalibration_",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0

drout <- alldrout[alldrout$idr==taskid,]
tolerance <- 1.5
drout <- drout[drout[,"rrinc"]/drout[,"inc"] > 1/tolerance*drout[,"targetdr"] & drout[,"rrinc"]/drout[,"inc"] < tolerance*drout[,"targetdr"], ] 


header <- c("inew", "ids","idr","targetprev","targetcoprev", "targetdr", "targetpt","DST", "rDSTall", names(unlist(mergedvalues)))
header <- append(header, paste0( rep(rtallynames, times=3*3*11), rep(0:10, each=length(rtallynames)),
                         "companion", rep(c("minimal","intermediate","optimal"), each=3*11*length(rtallynames)),
                         "barrier", rep(c("minimal","intermediate","optimal"), each=11*length(rtallynames)) ) )
                                                                 
if(!file.exists(paste0(location,"Resistance","_", targetpt,DST,"_",tasktag,".csv"))) { write(header,  file=paste0(location,"Resistance","_", targetpt,DST,"_",tasktag,".csv"), sep=",", ncol=length(header)) }
  
for (inew in 1:nrow(drout))
{
  iter <- unlist(c(inew,unlist(drout[inew,c("ids", "idr", "targetprev","targetcoprev","targetdr")]), targetpt, DST, rDSTall)) #will include these labels as part of returned output
  
  valuevect <- drout[inew, 5+(1:length(unlist(genericvalues)))] ; names(valuevect) <- names(unlist(genericvalues))
  v <- 0; for (pname in names(genericvalues))
  { genericvalues[[pname]] <- unlist(valuevect[v+(1:length(genericvalues[[pname]]))]); #names(genericvalues[[pname]]) <- names()
    v <- v + length(genericvalues[[pname]]) # move forward to start of next par vector in sampled values
  }    
  
  state <- drout[inew,drsetup$statenames]
  newstate <- numeric(length(novelsetup$statenames)); names(newstate) <- novelsetup$statenames
  newstate[drsetup$statenames] <- unlist(state)
  
  iresult <- numeric(0) #will become full row of results for this dr row (inew), for all elements and levels

  for (companion in c("minimal","intermediate","optimal")) for (barrier in c("minimal","intermediate","optimal"))
  {
    print(paste0("Evaluating TRP variation: companion ",companion,", barrier ",barrier," for Simulation #", inew," of ",nrow(drout)," for ",targetpt,DST,tasktag))
    valueset <- sampleTRP(mergedvalues = genericvalues, targetpt = targetpt, DST = DST, 
                          minimals=c("companion","barrier")[which(c(companion,barrier)=="minimal")], 
                          optimals=c("companion","barrier")[which(c(companion,barrier)=="optimal")]) 
    s_cr <- valueset$cres[1]; r_cr <- valueset$cres[2]
    novelstate <- newstate
    for (name in novelsetup$statenames) if (length(grep("^S", name))==0 & length(grep("^C", name))==0 )
    {  
      if (length(grep("+c+",name))==1 )
      {
        if (length(grep("+.Rr+", name))==1 ) 
        {    novelstate[name] <- r_cr * newstate[str_replace(string = name, pattern = "c", replacement = "")]
        } else novelstate[name] <- s_cr * max(newstate[str_replace(string = name, pattern = "c", replacement = "")], newstate[str_replace(string = name, pattern = "c", replacement = "0")], na.rm = TRUE)
      } else 
      {
        if (length(grep("+.Rr+", name))==1 )
        { novelstate[name] <- (1-r_cr) * newstate[name]
        } else novelstate[name] <- (1-s_cr) * newstate[name]
      }
    }
      
    parset <- create.pars(setup = novelsetup, values = valueset, T, T, T)
      
    outset <- ode(y=unlist(novelstate), times=0:10, func=dxdt, parms=parset$fullpars, do.tally=TRUE, method="adams")
      
    iresult <- append(iresult, c(t(outset[,rtallynames])))
  }
  write(unlist(c(iter, valuevect, iresult)), file=paste0(location,"Resistance","_", targetpt,DST,"_",tasktag,".csv"), sep=",", append=TRUE, ncol=length(header))
}
