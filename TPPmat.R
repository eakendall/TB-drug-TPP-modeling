library('deSolve')

# set up compartments
setup.model <- function()
{
  regimens <- c("s", "r", "n"); N <- 7; lengths <- c(2,1,1,2,3,3,6)/12
  Tnames <- c("S", "Ln", "An", "Ti", paste0("T", rep(regimens, times=N), rep(1:N, each=length(regimens))), "R", "C", "Lp", "Ap")
  Rnames <- c("R0", "Rc", "Rn", "Rcn", "Rr", "Rrc", "Rrn", "Rrcn") 
  Hnames <- c("Hn", "Hp")
  statenames <- paste0(rep(Tnames, times=length(Rnames)*length(Hnames)), ".", rep(Rnames, each=length(Tnames), times=length(Hnames)), ".", rep(Hnames, each=length(Tnames)*length(Rnames)))
  setup = list("regimens"=regimens, "Tperiods"=N, "Tperiod.lengths"=lengths, "Tnames"=Tnames, "Rnames"=Rnames, "Hnames"=Hnames, "statenames"=statenames)
  return(setup)
}


# define fixed parameters
## WILL WANT TO MAKE A SEPARATE FUNCTION FOR THE PARS THAT VARY IN THE NOVEL REGIMEN'S TPP, BUT FOR NOW EVERYTHING IS HERE IN ONE PLACE
setup.pars <- function(setup)
{
  if (missing(setup)) setup <- setup.model()
  p <- list()
  p$hivrate <- 0.001
  p$reactrate <- c(0.001,0.1); p$rapidprog <- c(0.1,0.5); p$mort <- c(0.03,0.1); p$tbmort <- c(0.2,0.4); 
      names(p$reactrate) <- names(p$rapidprog) <- names(p$mort) <- names(p$tbmort) <- setup$Hnames  #by HIV status, made up numbers for now
  p$selfcurerate <- 0.2; p$relapserate <- 1; p$latreduct <- 0.5
  p$dxrate <- c(1,2); p$DSTrif <- c(0.1,1); names(p$dxrate) <- names(p$DSTrif) <- c("An","Ap") #by new or previously treated
  p$availability <- 1 # MAY NEED TO MAKE THIS TIME-DEPENDENT FOR SCALE-UP
  p$DSTnew <- c(0,1); names(p$DSTnew) <- c("c","n") #companion drugs, novel drug (doesn't depend on rif or retreatment status)
  p$eligibility <- c(1,0.9); names(p$eligibility) <- setup$Hnames #by HIV status
  p$ltfurate <- c(0.01, 0.01, 0.01); names(p$ltfurate) <- setup$regimens
  p$targetpop <- c(1,1); names(p$targetpop) <- c("DS", "MDR") #will change for MDR=(0,1), DS=(1,0), or panTB=(1,1)
  
  fail_s <- c(0.02,0.5); fail_s <- rep(fail_s, each=4) #fraction ineffective on SOC DS regimen, by resistance
  fail_r <- 0.12; fail_r <- rep(fail_r, 8) #fraction ineffective on SOC MDR regimen
  fail_n <- c(fail_s[1:3] + c(0, 0.08, 0.5), 1); fail_n <- rep(fail_n, times=2)
  p$failmat <- rbind(fail_s, fail_r, fail_n); dimnames(p$failmat)= list(setup$regimens, setup$Rnames)
  
  p$transmissibility <- c(1,0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8); names(p$transmissibility) <- setup$Rnames #need to decide details of fitness costs and sampling
                        
  p$acqresmat <- array(0,dim=c(8,8,3)); dimnames(p$acqresmat)=list(setup$Rnames, setup$Rnames, setup$regimens) # from old resistance (down) to new resistance (across), by regimen 
  acqres_s <- 0.005; p$acqresmat[,,"s"] <- rbind(cbind(0,0,0,0, acqres_s * diag(4)), array(0, dim=c(4,8)))  
  p$acqresmat[,,"r"] <- array(0,dim=c(8,8)) # currently no acqres with regimen r
  acqres_n <- t(array(c( 0, 0.1, 0.05, 0.005,
                        0, 0, 0.4, 0.04, 
                        0, 0, 0, 0.1, 
                        0, 0, 0, 0), dim=c(4,4))); dimnames(acqres_n) <- list("R1"=c("0","c","n","cn"), "R2"=c("0","c","n","cn"))
  p$acqresmat[,,"n"] <- rbind(cbind(acqres_n,0,0,0,0), cbind(0,0,0,0, acqres_n))
  
  p$relapse_optimal <- c(0.04,0.12,0.02); names(p$relapse_optimal) <- setup$regimens #for fully susceptible and treatment completed
  p$relapse_resistance <- array(c( rep(p$relapse_optimal['s'],4), rep(0.6, 4), #translates 50% failure, 30% relapse, 20% cure into probabily of non-failures who relaspe
                                 rep(p$relapse_optimal['r'],8),  #again assuming resistance doesn't affect outcomes of SOC MDR regimen
                                 rep(c(p$relapse_optimal['n'], p$relapse_optimal['n']+0.06, 0.6, 1 ), 2) ), dim=c(8,3)) # adds ~6% relapse to total, or ~7% after failures and acqres removed
  dimnames(p$relapse_resistance) <- list(setup$Rnames, setup$regimens)
  
  p$durations <- c(6, 18, 4)/12; names(p$durations) <- setup$regimens #by regimen
  p$relapse246 <- c(0.3, 0.12, 0); names(p$relapse246) <- c("2mo","4mo","6mo") #extra relapse added by thirds of course completed (will interpolate between these)
  
  p <- c(p,setup)
  
  return(p)
}


# make matrix of non-dynamic state transitions (transitions are from first to second state, and diagonal is used for mortality)
makemat <- function(pars)
{
  if(missing(pars)) { pars <- setup.pars() }
  with(pars, 
  {mat <- array(0, dim=c(8+Tperiods*length(regimens), 8+Tperiods*length(regimens), 8, 8, 2, 2)); dimnames(mat) <- list("1T"=Tnames, "2T"=Tnames, "1R"=Rnames, "2R"=Rnames, "1H"=Hnames, "2H"=Hnames)
  
  ## HIV infection 
  for (jr in Rnames) for (jt in Tnames) {  mat[jt, jt, jr, jr,"Hn","Hp"] <- mat[jt, jt, jr, jr,"Hn","Hp"] + hivrate }
   
  # TB reactivation, mortality (background and TB-related, placed as ins on diagonal for now), self cure, and relapse:
  for (jr in Rnames) for (jh in Hnames) {
    mat["Ln","An",jr, jr, jh, jh] <- mat["Ln","An",jr, jr, jh, jh] + reactrate[jh]; 
    mat["Lp","Ap",jr, jr, jh, jh] <- mat["Lp","Ap",jr, jr, jh, jh] + reactrate[jh]  #reactivation, at hiv-dependent rate:
    diag(mat[,,jr,jr,jh,jh]) <- diag(mat[,,jr,jr,jh,jh]) + mort[jh]
    diag(mat[c("An","Ap","Ti"),c("An","Ap","Ti"),jr,jr,jh,jh]) <- diag(mat[c("An","Ap","Ti"),c("An","Ap","Ti"),jr,jr,jh,jh]) + tbmort[jh]
    for (jt in c("An","Ap","Ti")) { mat[jt,"C",jr,"R0",jh,jh] <- mat[jt,"C",jr,"R0",jh,jh] + selfcurerate } # note all cures go back to R0 to simplify later computations
    mat["R","Ap",jr,jr,jh,jh] <- mat["R","Ap",jr,jr,jh,jh] + relapserate
  }
  
  ## TB diagnosis and treatment
  for (jh in 1:2)
  {
    # tb is diagnosed and rif DST is performed depending on An vs Ap, then regimen is chosen according to MDR diagnosis, new regimen's target population, new regimen DST use, and medical eligibility by HIV status
    startmat <- array(0,dim=c(2,8,3)); dimnames(startmat) <- list("thist"=c("An","Ap"), Rnames, regimens) # probability of starting new regimen, by n/p treamtent history (down) and resistance (across) 
    startmat[,,"n"] <- cbind( 
    outer( targetpop[1] * c(1,1) * availability * eligibility[jh] , c(1, 1-DSTnew[1], 1-DSTnew[2], 1-max(DSTnew)) ) , # fraction of non-MDR (Rnames 1:4) that get novel regimen
    outer( (targetpop[1]*(1-DSTrif) + targetpop[2]*DSTrif) * availability * eligibility[jh] , c(1, 1-DSTnew[1], 1-DSTnew[2], 1-max(DSTnew)) ) # fraction of undiagnosed MDR can get new DS regimen, or diagnosed MDR can get new MDR regimen
    ) 
    startmat[,1:4,"s"] <- 1 - startmat[,1:4,"n"]; startmat[,5:8,"s"] <- (1-DSTrif) * (1 - startmat[,5:8,"n"]) 
    startmat[,5:8,"r"] <- DSTrif * (1 - startmat[,5:8,"n"])
    
    # now assign to ineffective treatment, or pending relapse with acquired resistance, or on first phase of effective treatment:
    for (it in c("An", "Ap")) 
    {
      for (nreg in 1:3)
      {
        mat[it, "R", , , jh, jh]  <- mat[it, "R", , , jh, jh] + dxrate[it] * startmat[it, ,nreg] * acqresmat[,,nreg] #acquired resistance moves to pending relapse for new strain
        diag( mat[it, "Ti", , , jh, jh]) <- diag( mat[it, "Ti", , , jh, jh]) + dxrate[it] * startmat[it, ,nreg] * failmat[nreg, ] #failures go to ineffective treatment (with ongoing infectiousness and increased mortality risk)
        diag( mat[it, c("Ts1","Tr1","Tn1")[nreg], , , jh, jh] ) <- diag( mat[it, c("Ts1","Tr1","Tn1")[nreg], , , jh, jh] ) + (dxrate[it] * startmat[it, , nreg]*(1 - failmat[nreg, ]- rowSums(acqresmat[,,nreg]))) #for the remainder, treatment is initially effective (i.e. outcome will be cure vs relapse if they complete at least 2 months)
      }
    }
        
    # determine whether  relapse or cure depending on effectiveness of regimen and fraction completed. Assuming losses occur in the middle of each time period on average.
    relapsefracs <- function(period) # use relapse %s at 2, 4, and 6 months (defined in pars) for a 6 month regimen, and interpolate linearly 2--4 and 4--6 to get relapse % after any fraction of treatment course
    {
        fraction_completed <- (sum(Tperiod.lengths[1:period-1])+1/2*Tperiod.lengths[period])/durations
        relapse <- array(0, dim=c(8,3)); dimnames(relapse) <- list(Rnames, regimens)
        relapse[,fraction_completed < 1/3] <- 1 
        relapse[,fraction_completed >= 1/3 & fraction_completed < 2/3] <- t(t(relapse_resistance + relapse246[2]) + (relapse246[1]-relapse246[2])*(2/3-fraction_completed))[, fraction_completed >= 1/3 & fraction_completed < 2/3]
        relapse[,fraction_completed >= 2/3 & fraction_completed < 1] <- t(t(relapse_resistance + relapse246[3]) + (relapse246[3]-relapse246[3])*(1-fraction_completed))[, fraction_completed >= 2/3 & fraction_completed < 1]
        relapse[,fraction_completed >= 1] <- relapse_resistance[, fraction_completed >= 1]
        relapse[relapse>1] <- 1 #if sum exceeds 100%, set to 100% relapse
        return(relapse) 
    }
    
    # once on effective treatment, progress through treatment to either relapse or cure (or to active if lost to follow up during first 2 months), including a monthly rate of loss to follow up 
    for (jr in Rnames) for (reg in regimens)
    {
      periods <- paste0("T", reg, 1:Tperiods)
      for (jt in 1:Tperiods)
      {
        relapses <- relapsefracs(period=jt)
        if (jt==1) 
        { mat[periods[1],c(periods[2],"Ap"),jr,jr,jh,jh] <- mat[periods[1],c(periods[2],"Ap"),jr,jr,jh,jh] + 1/Tperiod.lengths[1]*(c(0,1) + c(1,-1)*((1-ltfurate[reg])^Tperiod.lengths[1])) 
        } else 
          if (sum(Tperiod.lengths[1:jt]) < durations[reg]) 
        { mat[periods[jt],c(periods[jt+1],"C","R"),jr,jr,jh,jh] <- 
                mat[periods[jt],c(periods[jt+1],"C","R"),jr,jr,jh,jh] + 1/Tperiod.lengths[jt]*c( (1-ltfurate[reg])^Tperiod.lengths[jt], (1-(1-ltfurate[reg])^Tperiod.lengths[jt])*c(1-relapses[jr, reg], relapses[jr, reg])) 
        } else 
          if (sum(Tperiod.lengths[1:jt+1]) <= durations[reg] | jt==Tperiods) 
        { mat[periods[jt],c("C","R"),jr,jr,jh,jh] <- 
                mat[periods[jt],c("C","R"),jr,jr,jh,jh] + 1/Tperiod.lengths[jt]*c(1-relapses[jr, reg], relapses[jr, reg]) 
        }
      }
    }
  }
  return(mat)
  })
}




# dynamic model transitions

dxdt <- function(t, state, params, do.tally=FALSE)
{
 with(params, {
  if (length(mat) != length(state)^2) {stop("Error: Initial-state and transition-matrix size mismatch.")}
  
  statemat <- array(state, dim=c(length(Tnames), length(Rnames), length(Hnames))); dimnames(statemat) <- list(Tnames, Rnames, Hnames)
  FOI <- apply(statemat[c("An","Ap","Ti"),,], 2, sum) # FOI by strain
  
  # infection
  if (sum(statemat[c("S","C"), 2:8, ]) >0 ) { stop("Error: Some susceptibles have drug resistance and won't be included in infection events.") }
  for (jh in Hnames)
  {
    mat["S","Ln","R0",,jh,jh] <- mat["S","Ln","R0",,jh,jh] + FOI*transmissibility*(1-rapidprog[jh])
    mat["S","An","R0",,jh,jh] <- mat["S","An","R0",,jh,jh] + FOI*transmissibility*rapidprog[jh]
    mat["C","Lp","R0",,jh,jh] <- mat["C","Lp","R0",,jh,jh] + FOI*transmissibility*(1-rapidprog[jh])
    mat["C","Ap","R0",,jh,jh] <- mat["C","Ap","R0",,jh,jh] + FOI*transmissibility*rapidprog[jh]
  
    # superinfection
    for (jr in Rnames) 
    {
      mat["Ln","Ln",,jr,jh,jh] <- mat["Ln","Ln",,jr,jh,jh] + FOI*transmissibility[jr]*(1-rapidprog[jh]*latreduct)*transmissibility[jr]/(transmissibility + transmissibility[jr])
      mat["Ln","An",,jr,jh,jh] <- mat["Ln","An",,jr,jh,jh] + FOI*transmissibility[jr]*rapidprog[jh]*latreduct
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
    tally[grep("^A.[.]R[0cn]",statenames),"rrinc"] <-  # from active directly to relapse (when treatment starts) happens ony when resistance is acquired, so is a good place to count new rif resistance
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
    tally[grep("^A.[.]Rr",statenames),"nDSTs"] <- # if rif R 
      apply( squaremat[ grep("^A.[.]Rr",statenames), c(grep("^T",statenames),grep("^R",statenames)) ] * 
               ( (1-DSTrif)*targetpop[1] + DSTrif*targetpop[2] ) * availability * max(DSTnew) * rep(eligibility, each=length(grep("^A.[.]Rr",statenames))/2) , 1, sum )
    
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



#for debugging:
#dxdt(1, state, params, do.tally=FALSE)

#advance forward a bit:
state <- ode(state, seq(0,10,by=0.1), dxdt, params)[101,2:465]


runtoequilib(init, pars, ...)

addmdr( )


params <- pars(); params$mat <- make.mat(params)
