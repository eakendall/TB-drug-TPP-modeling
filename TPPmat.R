regimens <- c("s", "r", "n"); Tperiods <- 7; 
Tnames <- c("S", "Ln", "An", "Ti", paste0("T", rep(regimens, times=Tperiods), rep(1:Tperiods, each=length(regimens))), "R", "C", "Lp", "Ap")
Rnames <- c("R0", "Rc", "Rn", "Rcn", "Rr", "Rrc", "Rrn", "Rrcn") 
Hnames <- c("Hn", "Hp")

mat <- array(0, dim=c(8+Tperiods*length(regimens), 8+Tperiods*length(regimens), 8, 8, 2, 2)); dimnames(mat) <- list("1T"=Tnames, "2T"=Tnames, "1R"=Rnames, "2R"=Rnames, "1H"=Hnames, "2H"=Hnames)


hivrate <- 0.001
reactrate <- c(0.001,0.1); rapidprog <- c(0.1,0.5); mort <- c(0.03,0.1); tbmort <- c(0.2,0.4); 
    names(reactrate) <- names(rapidprog) <- names(mort) <- names(tbmort) <- Hnames  #by HIV status, made up numbers for now
selfcurerate <- 0.2; relapserate <- 1
dxrate <- c(1,2); DSTrif <- c(0.1,1); names(dxrate) <- names(DSTrif) <- c("An","Ap") #by new or previously treated
availability <- 1 # MAY NEED TO MAKE THIS TIME-DEPENDENT FOR SCALE-UP
DSTnew <- c(0,1) #companion drugs, novel drug (doesn't depend on rif or retreatment status)
eligibility <- c(1,0.9) #by HIV status
targetpop <- c(1,1) #will change for MDR=(0,1), DS=(1,0), or panTB=(1,1)
fail_s <- c(0.02,0.5); fail_s <- rep(fail_s, each=4) #fraction ineffective on SOC DS regimen, by rif resistance
fail_r <- 0.12; fail_r <- rep(fail_r, 8) #fraction ineffective on SOC MDR regimen
fail_n <- c(fail_s[1:3] + c(0, 0.07, 0.5), 1); fail_n <- rep(fail_n, times=2)
acqres_s <- 0.005 #fraction acquiring rif resistance on SOC DS regimen
# currently no acqres with regimen r
acqres_n <- t(array(c( 0, 0.1, 0.05, 0.005,
                        0, 0, 0.4, 0.04, 
                        0, 0, 0, 0.1, 
                        0, 0, 0, 0), dim=c(4,4)))


tally <- 0
# events to tally -- not yet implemented
if (tally==1) 
{
  tbcases <- 0
  tbdeaths <- 0
  rxtime <- c(0,0,0)
  riftests <- 0
  newtests <- 0
}
  

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
  for (jt in c("An","Ap","Ti")) { mat[jt,"C",jr,jr,jh,jh]) <- mat[jt,"C",jr,jr,jh,jh]) + selfcurerate }
  mat["R","Ap",jr,jr,jh,jh]) <- mat["R","Ap",jr,jr,jh,jh]) + relapserate
}


## Diagnosis and treatment initiation
for (jh in 1:2)
{
  # tb is diagnosed and rif DST is performed depending on An vs Ap, regimen is chosen according to MDR diagnosis, new regimen's target population, new regimen DST use, and medical eligibility by HIV status
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
  
  
  # some transition to failure as opposed to initial treatment response: (will redefine cures and nonres relapses %s to be _after_ failures as well as acq res are removed)
  failmat <- rbind(fail_s, fail_r, fail_n); dimnames(failmat)= list(regimens, Rnames)
  
  # then create an array for the probabilities of acquiring resistance when treatment starts (and going immediately to pending relapse for the new strain):
    # from old resistance (down) to new resistance (across), by regimen 
  acqresmat <- array(0,dim=c(8,8,3)); dimnames(acqresmat)=list(Rnames, Rnames, regimens)
  acqresmat[,,"n"] <- rbind(cbind(acqres_n,0,0,0,0), cbind(0,0,0,0, acqres_n))
  acqresmat[,,"s"] <- rbind(cbind(0,0,0,0, acqres_s * diag(4)), array(0, dim=c(4,8)))  
  acqresmat[,,"r"] <- array(0,dim=c(8,8))

  # now assign to ineffective treatment, or pending relapse with acquired resistance, or on first phase of effective treatment:
  for (it in c("An", "Ap")) 
  {
    for (reg in 1:3)
    {
      mat[it, "R", , , jh, jh]  <- mat[it, "R", , , jh, jh] + dxrate[it] * startmat[it, ,reg] * acqresmat[,,reg] #acquired resistance moves to pending relapse for new strain
      diag( mat[it, "Ti", , , jh, jh]) <- diag( mat[it, "Ti", , , jh, jh]) + dxrate[it] * startmat[it, ,reg] * failmat[reg, ] #failures go to ineffective treatment (with ongoing infectiousness and increased mortality risk)
      diag( mat[it, c("Ts1","Tr1","Tn1")[reg], , , jh, jh] ) <- diag( mat[it, c("Ts1","Tr1","Tn1")[reg], , , jh, jh] ) + (dxrate[it] * startmat[it, , reg]*(1 - failmat[reg, ]- rowSums(acqresmat[,,reg]))) #for the remainder, treatment is initially effective (i.e. outcome will be cure vs relapse if they complete at least 2 months)
    }
  }
      
  
## treatment course, with risk of ltfu after each time period, and with relapse vs cure depending on fraction of course completed



  
}





# births/ new adults
# born HIV negative



# time varying: 
infection