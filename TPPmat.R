library('deSolve')

# define fixed parameters; called by create.pars and carried forward in other functions
## !!  WILL WANT TO MAKE A SEPARATE FUNCTION FOR THE PARS THAT VARY IN THE NOVEL REGIMEN'S TPP, BUT FOR NOW EVERYTHING IS HERE IN ONE PLACE 

set.values <- function()
{
  values <- list()
  values <- within(values, {
  selfcurerate <- 0.2; relapserate <- 1; latreduct <- 0.5
  hivrate <- 0.001
  reactrate <- c(0.001,0.1); rapidprog <- c(0.1,0.5); mort <- c(0.03,0.1); tbmort <- c(0.2,0.4); #by HIV status, made up numbers for now
  beta <- 8; transmissibility <- 0.8
  
  months_s <- 6; months_r <- 18
  fail_s <- c(0.02,0.5); fail_r <- 0.12 # assumes second line outcomes are independent of resistance (and first line outcomes are independent of companion drug resistance)
  acqres_s <- 0.005; 
  relapse_s <- 0.04; relapse_r <- 0.12; mdrrelapse_s <- 0.6 #mdr relapse is of the 50% of rif resistance that doesn't fail outright
  dxrate <- c(0.5,1); DSTrif <- c(0.1,1); names(dxrate) <- names(DSTrif) <- c("An","Ap") 
  ltfurate_s <- 0.01; ltfurate_r <- 0.01
  relapse246 <- c(7.5, 3, 1); names(relapse246) <- c("2mo","4mo","6mo") #extra relapse added by thirds of course completed (will interpolate between these) -- multiplicative by regimen efficacy
  
  targetpop <- c(1,1); names(targetpop) <- c("DS", "MDR") #will change for MDR=(0,1), DS=(1,0), or panTB=(1,1)
  DSTnew <- c(0,1); names(DSTnew) <- c("c","n") #companion drugs, novel drug (doesn't depend on rif or retreatment status)
  months_n <- 4
  eligibility <- c(1,0.9)  #by HIV status
  availability <- 1 # MAY NEED TO MAKE THIS TIME-DEPENDENT FOR SCALE-UP !!
  ltfurate_n <- 0.01
  extrafail_n <- c(0, 0.08, 0.5) #for res to (-, c, n)
  acqres_n <- t(array(c( 0, 0.1, 0.05, 0.005, # down is starting resistance (-, c, n, cn), across is acquired pattern (-, c, n, cn) after novel regimen treatment
                         0, 0, 0.4, 0.04, 
                         0, 0, 0, 0.1, 
                         0, 0, 0, 0), dim=c(4,4)))
  relapse_n <- 0.02; resrelapse_n <- c(relapse_n+0.06, 0.6, 1 ) #risk of relapse (absolute, not additional) with novel regimen if resistance to c, n, and cn !!
  
  })
  
  return(values)
}


# sets up compartments; called by create.pars
setup.model <- function(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
{
  if ((treatSL & !DRera) | (treatnovel & !treatSL) | (treatnovel & !DRera)) stop("error: unexpected combination of available regimens")
  regimens <- c("s","r","n")
  if (treatnovel) { Rnames <- c("R0", "Rc", "Rn", "Rcn", "Rr", "Rrc", "Rrn", "Rrcn"); usereg <- c(1,2,3) } else
    if (treatSL) { Rnames <- c("R0", "Rr"); usereg <- c(1,2) } else 
      if (DRera) { Rnames <- c("R0", "Rr"); usereg <- 1 } else 
      { Rnames <- "R0"; usereg <- 1 }
  N <- 7; lengths <- c(2,1,1,2,3,3,6)/12
  Tnames <- c("S", "Ln", "An", "Ti", paste0("T", rep(regimens[usereg], times=N), rep(1:N, each=length(usereg))), "R", "C", "Lp", "Ap")
  
  Hnames <- c("Hn", "Hp")
  statenames <- paste0(rep(Tnames, times=length(Rnames)*length(Hnames)), ".", rep(Rnames, each=length(Tnames), times=length(Hnames)), ".", rep(Hnames, each=length(Tnames)*length(Rnames)))
  setup = list("DRera"=DRera, "treatSL"=treatSL, "treatnovel"=treatnovel, "regimens"=regimens, "usereg"=usereg, "Tperiods"=N, "Tperiod.lengths"=lengths, "Tnames"=Tnames, "Rnames"=Rnames, "Hnames"=Hnames, "statenames"=statenames)
  return(setup)
}


# now combine values and model setup into properly-formatted parameters
# returns list of lists: parameters (pars, for use) and stand-alone values (in case need to change some values or setup later)
create.pars <- function(setup, values, DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
{
  if (missing(setup)) setup <- setup.model(DRera, treatSL, treatnovel)
  if (missing(values)) values <- set.values() 
  pars <- values
  
  with(setup, {
    pars <- within(pars, {
      
      names(reactrate) <- names(rapidprog) <- names(mort) <- names(tbmort) <- Hnames  
      
      transmissibility <- c(1,rep(transmissibility, length(Rnames)-1)); names(transmissibility) <- Rnames #need to decide details of fitness costs and sampling
      
      # then adjust for whether regimens are in use (i.e.if not using FL, no diagnosis and treatment; if not using SL, no rif diagnosis (and therefore no second-line treatment); 
      # and if not using novel regimen (which can be used in absence of novel DST in the model), availability of novel regimen is 0.
      if (DRera==FALSE) acqres_s <- 0; if (treatSL==FALSE) DSTrif <- c(0,0); if (treatnovel==FALSE) availability <- 0 
      
      names(eligibility) <- Hnames
      ltfurate <- c(ltfurate_s, ltfurate_r, ltfurate_n); names(ltfurate) <- regimens
      
      #translate above parameters into outcome arrays (shouldn't need to edit this part unless model structure changes):
      fail_s <- c(rep(fail_s[1], max(1,4*treatnovel)), rep(fail_s[2], DRera*max(1,4*treatnovel)))
      fail_n <- c(fail_s[1] + extrafail_n, 1) #for res to (-, c, n), cn
      failmat <- array(rbind(fail_s, fail_r, fail_n)[,1:length(Rnames)], dim=c(3,length(Rnames))); dimnames(failmat) <- list(regimens, Rnames) 
      
      acqresmat <- array(0,dim=c(length(Rnames),length(Rnames),length(regimens))); dimnames(acqresmat)=list(Rnames, Rnames, regimens) # from old resistance (down) to new resistance (across), by regimen 
      acqresmat[grep("^R[0cn]",Rnames), grep("^Rr",Rnames),"s"] <-  acqres_s * diag(max(1, 4*treatnovel))
      acqresmat[,,"r"] <- array(0,dim=c(length(Rnames),length(Rnames))) # currently no acqres with regimen r
      if (treatnovel) { acqresmat[,,"n"] <- rbind(cbind(acqres_n,0,0,0,0), cbind(0,0,0,0, acqres_n))[1:length(Rnames), 1:length(Rnames)] }
      
      relapse_optimal <- c(relapse_s,relapse_r,relapse_n); names(relapse_optimal) <- regimens #for fully susceptible and treatment completed
      
      relapse_resistance <- array(c( rep(relapse_s,length(grep("^R[0cn]",Rnames))), rep(mdrrelapse_s, length(grep("^Rr",Rnames))), #translates 50% failure, 30% relapse, 20% cure into probabily of non-failures who relaspe
                                     rep(relapse_r,length(Rnames)),  #again assuming resistance doesn't affect outcomes of SOC MDR regimen
                                     rep(c(relapse_n, resrelapse_n), 2)[1:length(Rnames)]), dim=c(length(Rnames),3) ) # adds ~6% relapse to total, or ~7% after failures and acqres removed
      dimnames(relapse_resistance) <- list(Rnames, regimens)
      
      durations <- c(months_s, months_r, months_n)/12; names(durations) <- regimens #by regimen
      })
    pars <- c(pars,setup)
    pars$mat <- makemat(pars=pars)
  
    return(list("fullpars"=pars, "values"=values))
  })
  
}


# creates matrix of non-dynamic state transitions to be added to parameters (transitions are from first to second state, and diagonal is used for mortality)
# Intended to be called within create.pars; argument pars is the fullpars output of create.pars prior to adding mat
makemat <- function(pars) 
{
  with(pars, {
    mat <- array(0, dim=c(length(Tnames), length(Tnames), length(Rnames), length(Rnames), length(Hnames), length(Hnames))); dimnames(mat) <- list("1T"=Tnames, "2T"=Tnames, "1R"=Rnames, "2R"=Rnames, "1H"=Hnames, "2H"=Hnames)
    
    ## HIV infection 
    for (jr in Rnames) for (jt in Tnames) {  mat[jt, jt, jr, jr,"Hn","Hp"] <- mat[jt, jt, jr, jr,"Hn","Hp"] + hivrate }
     
    # TB reactivation, mortality (background and TB-related, placed as ins on diagonal for now), self cure (all to R0 resistance), and relapse:
    for (jr in Rnames) for (jh in Hnames) {
      mat["Ln","An",jr, jr, jh, jh] <- mat["Ln","An",jr, jr, jh, jh] + reactrate[jh]; 
        mat["Lp","Ap",jr, jr, jh, jh] <- mat["Lp","Ap",jr, jr, jh, jh] + reactrate[jh]  #reactivation, at hiv-dependent rate:
      diag(mat[,,jr,jr,jh,jh]) <- diag(mat[,,jr,jr,jh,jh]) + mort[jh]
      diag(mat[c("An","Ap","Ti"),c("An","Ap","Ti"),jr,jr,jh,jh]) <- diag(mat[c("An","Ap","Ti"),c("An","Ap","Ti"),jr,jr,jh,jh]) + tbmort[jh]
      mat["An","S",jr,"R0",jh,jh] <- mat["An","S",jr,"R0",jh,jh] + selfcurerate;  # note all cures go back to R0 to simplify later computations
        for (jt in c("Ap","Ti")) { mat[jt,"C",jr,"R0",jh,jh] <- mat[jt,"C",jr,"R0",jh,jh] + selfcurerate } 
      mat["R","Ap",jr,jr,jh,jh] <- mat["R","Ap",jr,jr,jh,jh] + relapserate
    }
    
    ## TB diagnosis and treatment
    for (jh in 1:length(Hnames))
    {
      # tb is diagnosed and rif DST is performed depending on An vs Ap, 
      # then regimen is chosen according to MDR diagnosis, new regimen's target population, new regimen DST use, and medical eligibility by HIV status
      # considers all resistance including those not included in model setup
      startmat <- array(0,dim=c(2,8,3)); dimnames(startmat) <- list("thist"=c("An","Ap"), c("R0", "Rc", "Rn", "Rcn", "Rr", "Rrc", "Rrn", "Rrcn"), regimens) # probability of starting new regimen, by n/p treamtent history (down) and resistance (across) 
      
      startmat[,,"n"] <- cbind( 
      outer( targetpop[1] * c(1,1) * availability * eligibility[jh] , c(1, 1-DSTnew[1], 1-DSTnew[2], 1-max(DSTnew)) ) , # fraction of non-MDR (Rnames 1:4) that get novel regimen
      outer( (targetpop[1]*(1-DSTrif) + targetpop[2]*DSTrif) * availability * eligibility[jh] , c(1, 1-DSTnew[1], 1-DSTnew[2], 1-max(DSTnew)) ) # fraction of undiagnosed MDR can get new DS regimen, or diagnosed MDR can get new MDR regimen
      ) 
      startmat[,1:4,"s"] <- 1 - startmat[,1:4,"n"]; startmat[,5:8,"s"] <- (1-DSTrif) * (1 - startmat[,5:8,"n"]) 
      startmat[,5:8,"r"] <- DSTrif * (1 - startmat[,5:8,"n"])
    
      # now assign to ineffective treatment, or pending relapse with acquired resistance, or on first phase of effective treatment: !! full %s are removed for both fail and acqres before the rest get started on treatment (or, another way to think about it is that acqres come from both failures and treatment responsives, proportionately)
      for (it in c("An", "Ap")) 
      {
        for (nreg in usereg)
        {
          mat[it, "R", , , jh, jh]  <- 
            mat[it, "R", , , jh, jh] + dxrate[it] * startmat[it,Rnames,nreg] * acqresmat[,,nreg] #acquired resistance moves to pending relapse for *new* strain (acq)
          
          if (length(Rnames)>1)
          {
            diag(mat[it, "Ti", , , jh, jh]) <- 
              diag(mat[it, "Ti", , , jh, jh]) + dxrate[it] * startmat[it, Rnames ,nreg] * failmat[nreg, ] #failures go to ineffective treatment (with ongoing infectiousness and increased mortality risk)
            diag(mat[it, grep("T[srn]1",Tnames)[nreg], , , jh, jh]) <- 
              diag(mat[it, grep("T[srn]1",Tnames)[nreg], , , jh, jh]) + dxrate[it] * startmat[it, Rnames, nreg]*(1 - failmat[nreg, ]- rowSums(acqresmat[,,nreg])) #for the remainder, treatment is initially effective (i.e. outcome will be cure vs relapse if they complete at least 2 months)
          } else
          {
            mat[it, "Ti", 1,1, jh, jh] <- 
              mat[it, "Ti", 1,1, jh, jh] + dxrate[it] * startmat[it, Rnames ,nreg] * failmat[nreg, ] #failures go to ineffective treatment (with ongoing infectiousness and increased mortality risk)
            mat[it, grep("T[srn]1",Tnames)[nreg], , , jh, jh] <- 
              mat[it, grep("T[srn]1",Tnames)[nreg], , , jh, jh] + (dxrate[it] * startmat[it, Rnames, nreg]*(1 - failmat[nreg, ]- sum(acqresmat[,,nreg]))) #for the remainder, treatment is initially effective (i.e. outcome will be cure vs relapse if they complete at least 2 months)            
          }
        }
      }
          
        # determine whether  relapse or cure depending on effectiveness of regimen and fraction completed. Assuming losses occur in the middle of each time period on average.
        relapsefracs <- function(period) # use relapse %s at 2, 4, and 6 months (defined in pars) for a 6 month regimen, and interpolate linearly 2--4 and 4--6 to get relapse % after any fraction of treatment course
        {
            fraction_completed <- (sum(Tperiod.lengths[1:period-1])+1/2*Tperiod.lengths[period])/durations
            relapse <- array(0, dim=c(length(Rnames),3)); dimnames(relapse) <- list(Rnames, regimens)
            relapse[,fraction_completed < 1/3] <- 1 
            relapse[,fraction_completed >= 1/3 & fraction_completed < 2/3] <- 
              t(t(relapse_resistance) * (relapse246[2] + (relapse246[1]-relapse246[2])*(2/3-fraction_completed)))[, fraction_completed >= 1/3 & fraction_completed < 2/3]
            relapse[,fraction_completed >= 2/3 & fraction_completed < 1] <- 
              t(t(relapse_resistance) * (relapse246[3]) + (relapse246[3]-relapse246[3])*(1-fraction_completed))[, fraction_completed >= 2/3 & fraction_completed < 1]
            relapse[,fraction_completed >= 1] <- relapse_resistance[, fraction_completed >= 1]
            relapse[relapse>1] <- 1 #if sum exceeds 100%, set to 100% relapse
            return(relapse) 
        }
        
        # once on effective treatment, progress through treatment to either relapse or cure (or to active if lost to follow up during first 2 months), including a monthly rate of loss to follow up 
        if (length(usereg) > 0) for (jr in Rnames) for (reg in regimens[usereg])
        {
          periods <- paste0("T", reg, 1:Tperiods)
          jt <- 1
          while (sum(Tperiod.lengths[1:jt]) < durations[reg])
          {
            relapses <- relapsefracs(period=jt)
            mat[periods[jt],c(periods[jt+1],"R"),jr,jr,jh,jh] <- # move on; or discontinue and relapse 
              mat[periods[jt],c(periods[jt+1],"R"),jr,jr,jh,jh] + 1/Tperiod.lengths[jt]* c( (1-ltfurate[reg])^(12*Tperiod.lengths[jt]), (1-(1-ltfurate[reg])^(12*Tperiod.lengths[jt]))*relapses[jr, reg])
            mat[periods[jt],"C",jr,"R0",jh,jh] <- # discontinue and cured (goes back to no resistance)
              mat[periods[jt],"C",jr,"R0",jh,jh] + 1/Tperiod.lengths[jt]* (1-(1-ltfurate[reg])^(12*Tperiod.lengths[jt])) * (1-relapses[jr, reg])
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




# function for ode solver, including dynamic model transitions
# fullpars argument is the fullpars component of create.pars (i.e. list of params incl mat, but doesn't include list of values)
# output is ode output, along with outcomes (as defined within function) at each time point
dxdt <- function(t, state, fullpars, do.tally=FALSE)
{
  with(fullpars, {
  if (length(mat) != length(state)^2) {stop("Error: Initial-state and transition-matrix size mismatch.")}
  
  statemat <- array(state, dim=c(length(Tnames), length(Rnames), length(Hnames))); dimnames(statemat) <- list(Tnames, Rnames, Hnames)
  if (length(Rnames)==1) { FOI = sum(statemat[c("An","Ap","Ti"),,]) * transmissibility * beta / sum(statemat) 
    } else  { FOI <- apply(statemat[c("An","Ap","Ti"),,], 2, sum) * transmissibility * beta / sum(statemat) }# FOI by strain
  
  # infection
#   if (length(Rnames)>1 && sum(statemat[c("S","C"), -1, ]) >0 ) { stop("Error: Some susceptibles have drug resistance and won't be included in infection events.") }
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

  ########### transform mat
  newmat <- aperm(mat, c(1,3,5,2,4,6)) # this will translate into a 2d matrix more easily
  squaremat <- array(newmat, dim=c(length(Tnames)*length(Rnames)*length(Hnames), length(Tnames)*length(Rnames)*length(Hnames))) #2d state1 to state2 transmition matrix
  dimnames(squaremat) <- list(statenames, statenames)
  
  
  ########## tally how much each current state contributes to outcomes of interest (then will multiply by state and append to output)
  
  outcomes <- c("prev", "inc", "rrinc", "relapses", "tbdeaths", "rrdeaths", "dxs", "rDSTs", "nDSTs", "rxtime_s", "rxtime_r", "rxtime_n")
  tally <- array(0,dim=c(length(Tnames)*length(Rnames)*length(Hnames), length(outcomes))); dimnames(tally) <- list(statenames, outcomes)
  
  if (do.tally==TRUE)
  {
    tally[c(grep("^A", statenames), grep("^T", statenames)), "prev"] <- rep(1, length(c(grep("^A", statenames), grep("^T", statenames))))
    
    tally[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames)),"inc"] <- 
      apply(squaremat[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames)),  grep("^A",statenames)], 1, sum)
  
    tally[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames)),"rrinc"] <- 
      apply(squaremat[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames)),  grep("^A.[.]Rr",statenames)], 1, sum)
    tally[grep("^A.[.]R[0cn]",statenames),"rrinc"] <-  # from active directly to relapse (when treatment starts) happens ony when resistance is acquired, so is a good place to count new rif resistance (but will produce small nonzero incidence at time zero)
      apply(squaremat[grep("^A.[.]R[0cn]",statenames),  grep("^R.Rr",statenames)], 1, sum)
    
    tally[grep("^R",statenames),"relapses"] <- 
      apply(squaremat[grep("^R",statenames),  grep("^A",statenames)], 1, sum)
    
    tally[c(grep("^A", statenames), grep("^Ti", statenames)), "tbdeaths"] <- rep(tbmort, each=(length(c(grep("^A", statenames), grep("^Ti", statenames)))/2))
    
    tally[c(grep("^A.[.]Rr", statenames), grep("^Ti[.]Rr", statenames)), "rrdeaths"] <- rep(tbmort, each=(length(c(grep("^A.[.]Rr", statenames), grep("^Ti[.]Rr", statenames)))/2))
    
    tally[grep("^A", statenames),"dxs"] <- apply( squaremat[ grep("^A", statenames), c(grep("^T", statenames), grep("^R", statenames)) ], 1, sum)
    
    tally[grep("^A", statenames),"rDSTs"] <- 
      apply( squaremat[ grep("^A", statenames), c(grep("^T",statenames),grep("^R",statenames)) ] * DSTrif, 1, sum) #alternates between An and Ap, so alternate DST coverage
             
    tally[grep("^A.[.]R[0cn]",statenames),"nDSTs"] <- # if not rif R
      apply( squaremat[ grep("^A.[.]R[0cn]",statenames), c(grep("^T",statenames),grep("^R",statenames)) ] * 
               targetpop[1] * availability * max(DSTnew) * rep(eligibility, each=length(grep("^A.[.]R[0cn]",statenames))/2) , 1, sum )
    if (length(grep("^A.[.]Rr",statenames))>0)
    {
      tally[grep("^A.[.]Rr",statenames),"nDSTs"] <- # if rif R 
        apply( squaremat[ grep("^A.[.]Rr",statenames), c(grep("^T",statenames),grep("^R",statenames)) ] * 
                 ( (1-DSTrif)*targetpop[1] + DSTrif*targetpop[2] ) * availability * max(DSTnew) * rep(eligibility, each=length(grep("^A.[.]Rr",statenames))/2) , 1, sum )
    }
    
    tally[grep("^Ts", statenames), "rxtime_s"] <- rep(1, length(grep("^Ts", statenames))) #doesn't include ineffective (Ti) months
    tally[grep("^Tr", statenames), "rxtime_r"] <- rep(1, length(grep("^Tr", statenames))) #doesn't include ineffective (Ti) months
    tally[grep("^Tn", statenames), "rxtime_n"] <- rep(1, length(grep("^Tn", statenames))) #doesn't include ineffective (Ti) months
  }
  
  tallied <- t(tally) %*% state ; names(tallied) <- outcomes 
  
  ## combine state vector and change matrix to get dxdt
  dxdt <- numeric(length(state))
  deaths <- state*c(diag(squaremat)); diag(squaremat) <- 0 #store deaths elsewhere
  dxdt <- dxdt + t(squaremat) %*% state - rowSums(squaremat*state) - deaths #ins minus outs minus deaths
    
  # births/ new adults (# born HIV negative and TB susceptible, or make a FOI-dependent latent TB fraction?)
  dxdt[1] <- dxdt[1] + sum(deaths)  
  
  # return dxdt to ode, and also return tally for tracking purposes
  return(list(dxdt, tallied))
 })
}



#advance forward a bit:
## stateplus includes outcomes tracking variables
advance <- function(state, t0, dxdt, addedt=1, reportsteps=1, fullpars) 
{
  return(ode(state, seq(t0,t0+addedt,by=0.1), dxdt, fullpars, do.tally=TRUE)[(10*addedt+1 - reportsteps):(10*addedt+1),])
}

  
# run to equilibrium without drugs; will need to adjust Rnames and resistance-related parameters, and remake mat
equilib <- function(state, func=dxdt, pars, tol=0.1)
{
  if(missing(pars)) pars <- create.pars(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
  
  if(missing(state)) { state <- with(pars$fullpars,  c(90000, 9900, 100, rep(0,length(statenames)-3))); names(state) <- pars$fullpars$statenames }
  
  statex <- ode(state, seq(0,20,by=0.1), func, pars$fullpars, do.tally=TRUE)[200:201,]; totaltime <- 20
  
  log <- c(0, state, rep(0, ncol(statex)-1-length(state))); names(log) <- colnames(statex)
  log <- rbind(log, statex[2,])
  
  while (max (statex[2,-1]-statex[1,-1])>tol)
  {
    statex <- advance(state=statex[2,2:(length(state)+1)], t0=totaltime, addedt=10, func, fullpars=pars$fullpars)
    totaltime <- totaltime + 10
    log <- rbind(log, statex[2,])
  }
  return(list("log"=log, "pars"=pars))
}
  
plotlog <- function(log)
{
  par(mfrow=c(2,2), mar=c(2,4,1,1))
  plot(log[,1], log[,"inc"], type='l') 
  plot(log[,1], log[,"prev"]) 
  plot(log[,1], log[,"tbdeaths"]) 
  plot(log[,1], log[,"dxs"]) 
  
}

addmdr<- function(eqb, acqres=0.005, DSTrif = c(0,0), rtrans=0.7, addedt=10)
{
  if(missing(eqb)) eqb <- equilib()
  
  oldstatenames <- eqb$pars$fullpars$statenames
  oldstate <- with(eqb,log[nrow(log),2:(length(oldstatenames)+1)])
  
  oldvalues <- eqb$pars$values
  newvalues <- oldvalues; newvalues$acqres_s <- acqres; newvalues$DSTrif <- DSTrif; newvalues$transmissibility <- rtrans
  
  if (sum(DSTrif)>0) treatSL <- TRUE else treatSL <- FALSE
  newpars <- create.pars(values=newvalues, DRera=TRUE, treatSL=treatSL, treatnovel=FALSE)
  
  newstate <- numeric(length(newpars$fullpars$statenames)); names(newstate) <- newpars$fullpars$statenames
  newstate[oldstatenames] <- oldstate
  
  mdrode <- advance(newstate, t0=0, dxdt, addedt, reportsteps = 10*addedt, fullpars=newpars$fullpars) 
  
  par(mfrow=c(2,2), mar=c(2,4,1,1))
  plot(mdrode[,"time"], mdrode[,"inc"])
  plot(mdrode[,"time"], mdrode[,"rrinc"])
  plot(mdrode[,"time"], mdrode[,"rrdeaths"])
  plot(mdrode[,"time"], mdrode[,"rDSTs"])
  
  return(list("mdrode"=mdrode, "pars"=newpars))
}

