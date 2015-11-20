library('deSolve')

# set up compartments
setup.model <- function()
{
  regimens <- c("s", "r", "n"); N <- 7; lengths <- c(2,1,1,2,3,3,6)
  Tnames <- c("S", "Ln", "An", "Ti", paste0("T", rep(regimens, times=Tperiods), rep(1:Tperiods, each=length(regimens))), "R", "C", "Lp", "Ap")
  Rnames <- c("R0", "Rc", "Rn", "Rcn", "Rr", "Rrc", "Rrn", "Rrcn") 
  Hnames <- c("Hn", "Hp")
  setup = list("regimens"=regimens, "Tperiods"=N, "Tperiod.lengths"=lengths, "Tnames"=Tnames, "Rnames"=Rnames, "Hnames"=Hnames)
  return(setup)
}

# define parameters
## WILL WANT TO MAKE A SEPARATE FUNCTION FOR THE PARS THAT VARY IN THE NOVEL REGIMEN'S TPP, BUT FOR NOW EVERYTHING IS HERE IN ONE PLACE
pars <- function(setup)
{
  if (missing(setup)) setup <- setup.model()
  p <- list()
  p$hivrate <- 0.001
  p$reactrate <- c(0.001,0.1); p$rapidprog <- c(0.1,0.5); p$mort <- c(0.03,0.1); p$tbmort <- c(0.2,0.4); 
      names(p$reactrate) <- names(p$rapidprog) <- names(p$mort) <- names(p$tbmort) <- setup$Hnames  #by HIV status, made up numbers for now
  p$selfcurerate <- 0.2; p$relapserate <- 1
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
  
  p$durations <- c(6, 18, 4); names(p$durations) <- setup$regimens #in months, by regimen
  p$relapse246 <- c(0.3, 0.12, 0); names(p$relapse246) <- c("2mo","4mo","6mo") #extra relapse added by thirds of course completed (will interpolate between these)
  
  p <- c(p,setup)
  
  return(p)
}



# make matrix of non-dynamic state transitions (transitions are from first to second state, and diagonal is used for mortality)
makemat <- function(pars)
{
  with(pars, 
  {mat <- array(0, dim=c(8+Tperiods*length(regimens), 8+Tperiods*length(regimens), 8, 8, 2, 2)); dimnames(mat) <- list("1T"=Tnames, "2T"=Tnames, "1R"=Rnames, "2R"=Rnames, "1H"=Hnames, "2H"=Hnames)
  mat
  
  ## HIV infection 
  for (jr in Rnames) for (jt in Tnames) {  mat[jt, jt, jr, jr,"Hn","Hp"] <- mat[jt, jt, jr, jr,"Hn","Hp"] + hivrate }
   
  ## TB natural history 
  #infections (resulting in latent or active) inluding superinfections go in main model
  # reactivation, mortality (background and TB-related, placed as ins on diagonal for now), self cure, and relapse:
  for (jr in Rnames) for (jh in Hnames) {
    mat["Ln","An",jr, jr, jh, jh] <- mat["Ln","An",jr, jr, jh, jh] + reactrate[jh]; 
    mat["Lp","Ap",jr, jr, jh, jh] <- mat["Lp","Ap",jr, jr, jh, jh] + reactrate[jh]  #reactivation, at hiv-dependent rate:
    diag(mat[,,jr,jr,jh,jh]) <- diag(mat[,,jr,jr,jh,jh]) + mort[jh]
    diag(mat[c("An","Ap","Ti"),c("An","Ap","Ti"),jr,jr,jh,jh]) <- diag(mat[c("An","Ap","Ti"),c("An","Ap","Ti"),jr,jr,jh,jh]) + tbmort[jh]
    for (jt in c("An","Ap","Ti")) { mat[jt,"C",jr,jr,jh,jh] <- mat[jt,"C",jr,jr,jh,jh] + selfcurerate }
    mat["R","Ap",jr,jr,jh,jh] <- mat["R","Ap",jr,jr,jh,jh] + relapserate
  }
  
  ## Diagnosis and treatment
  for (jh in 1:2)
  {
    # tb is diagnosed and rif DST is performed depending on An vs Ap, then regimen is chosen according to MDR diagnosis, new regimen's target population, new regimen DST use, and medical eligibility by HIV status
    # first create arrays for the fraction who start each regimen once diagnosed, by treatment history and by resistance:
    startmat <- array(0,dim=c(2,8,3)); dimnames(startmat) <- list("thist"=c("An","Ap"), Rnames, regimens)
    startmat[,,"n"] <- cbind( # probability of starting new regimen, by n/p treamtent history (down) and resistance (across) 
    # fraction of non-MDR (Rnames 1:4) that get novel regimen, by n/p treatment history
    outer( targetpop[1] * c(1,1) * availability * eligibility[jh] , 
     c(1, 1-DSTnew[1], 1-DSTnew[2], 1-max(DSTnew)) ) #(by treatment hist, new reg resistance)
    ,
    # fraction of MDR (Rnames 5:8) that get novel regimen, by n/p treatment history
    outer(   (targetpop[1]*(1-DSTrif) + targetpop[2]*DSTrif) * availability * eligibility[jh] , #undiagnosed MDR can get new DS regimen, or diagnosed MDR can get new MDR regimen
     c(1, 1-DSTnew[1], 1-DSTnew[2], 1-max(DSTnew)) ) #(by treatment hist, new reg resistance)
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
        
    
    ## progression through effective treatment to either relapse or cure (or to active if lost to follow up during first 2 months)
    
    # risk of loss to follow up during each time period, and with relapse vs cure depending on fraction of course completed (assuming losses occur in the middle of each time period on average)
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

    
    # progression through treatment, including a monthly rate of loss to follow up
    for (jr in Rnames) for (reg in regimens)
    {
      periods <- paste0("T", reg, 1:Tperiods)
      for (jt in 1:Tperiods)
      {
        relapses <- relapsefracs(period=jt)
        if (jt==1) 
        { mat[periods[1],c(periods[2],"Ap"),jr,jr,jh,jh] <- mat[periods[1],c(periods[2],"Ap"),jr,jr,jh,jh] + (c(0,1) + c(1,-1)*((1-ltfurate[reg])^Tperiod.lengths[1])) 
        } else 
          if (sum(Tperiod.lengths[1:jt]) < durations[reg]) 
        { mat[periods[jt],c(periods[jt+1],"C","R"),jr,jr,jh,jh] <- 
                mat[periods[jt],c(periods[jt+1],"C","R"),jr,jr,jh,jh] + c( (1-ltfurate[reg])^Tperiod.lengths[jt], (1-(1-ltfurate[reg])^Tperiod.lengths[jt])*c(1-relapses[jr, reg], relapses[jr, reg])) 
        } else 
          if (sum(Tperiod.lengths[1:jt+1]) <= durations[reg] | jt==Tperiods) 
        { mat[periods[jt],c("C","R"),jr,jr,jh,jh] <- 
            mat[periods[jt],c("C","R"),jr,jr,jh,jh] + c(1-relapses[jr, reg], relapses[jr, reg]) 
        }
      }
    }
  
  }
  return(mat)
  })
}




# dynamic model transitions

dxdt <- function(t, state, params)
{
  # infection
  # superinfection
  
  
 
  # dxdt  
  
  # births/ new adults (# born HIV negative and TB susceptible, or make a FOI-dependent latent TB fraction?)
  # events to tally -- not yet implemented. move to dynamic model? 

  if (tally==1) 
  {
    tbcases <- 0
    tbdeaths <- 0
    rxtime <- c(0,0,0)
    riftests <- 0
    newtests <- 0
  }
  
}

run.model(dxdt, init, pars, ...)


params <- pars(); params$mat <- make.mat(params)
