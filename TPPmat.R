library('deSolve')

# define fixed parameters; called by create.pars and carried forward in other functions

# will provide a list of "values" inputs for the function that implements novel regimen: 
# each labeled by the targetpop, the DST use, the TRP element varied (with "none" as one option), and whether the varied TRP element is minimal or optimal
samplenovel <- function(values, target="DS", DST=FALSE)
{
  if (target == "DS") {values$targetpop <- c(1,0)} else {values$targetpop <- c(0,1)}
  if (DST) {values$DSTnew[1:2] <- c(1,1)} else {values$DSTnew[1:2] <- c(0,0)}
    
  #set all TPR elements to middle value
  if (target=="DS") {poor_n <- 0.03} else {poor_n <- 0.06}
  if (target=="DS") {values$months_n <- 4} else {values$months_n <- 9}
  if (target=="DS") {values$cres[1:2] <- rep(0.03, 2)} else {values$cres[1:2] <-rep(0.1, 2)}  #or can make the two vector elements different, if companion resistance is correlated with rif resistance
  if (target=="DS") {barrierbase <- 0.05} else {barrierbase <- 0.008} 
  values$eligibility <- 1 - rep(0.02, 2) # can vary by HIV status !!
  if (target=="DS") {values$ltfurate_n <- values$ltfurate_s - 0.015/values$months_n} else {values$ltfurate_n <- values$ltfurate_s - 0.03/values$months_n}
      
  values$poor_n <- c(poor_n + c(0, 0.13, 0.5), 1)
  
  values$acqres_n <- t(array(c( 0, 0.1, barrierbase, 0.1*barrierbase, # down is starting resistance (-, c, n, cn), across is acquired pattern (-, c, n, cn) after novel regimen treatment
                                0, 0, 0, 8*barrierbase, 
                                0, 0, 0, 0.1, 
                                0, 0, 0, 0), dim=c(4,4))); values$acqres_n[values$acqres_n>1] <- 1
  
  elementnames <- c("effectiveness", "duration", "companion", "barrier", "exclusions", "exclusions_hiv", "tolerability"); levelnames <- c("minimal", "optimal")
  
  levels <- as.list(levelnames); names(levels) <- levelnames
  for (name in levelnames) levels[[name]] <- values
  TRP <- as.list(elementnames); names(TRP) <- elementnames
  for (name in elementnames) TRP[[name]] <- levels
  TRP$none <- list(); TRP$none$none <- values
  # input minimal and optimal TRP ranges here
  if (target=="DS")
  {
    poor_n <- 0.06; TRP$effectiveness$minimal$poor_n <- c(poor_n + c(0, 0.13, 0.5), 1)
    poor_n <- 0; TRP$effectiveness$minimal$poor_n <- c(poor_n + c(0, 0.13, 0.5), 1)
    
    TRP$duration$minimal$months_n <- 6
    TRP$duration$optimal$months_n <- 3
    
    TRP$companion$minimal$cres <- rep(0.1,2)
    TRP$companion$optimal$cres <- rep(0,2)
    
    barrierbase <- 0.05; TRP$barrier$minimal$acqres_n <- t(array(c( 0, 0.1, barrierbase, 0.1*barrierbase, # down is starting resistance (-, c, n, cn), across is acquired pattern (-, c, n, cn) after novel regimen treatment
                                                                        0, 0, 0, 8*barrierbase, 
                                                                        0, 0, 0, 0.1, 
                                                                        0, 0, 0, 0), dim=c(4,4))); TRP$barrier$minimal$acqres_n[TRP$barrier$minimal$acqres_n>1] <- 1
    barrierbase <- 0; TRP$barrier$minimal$acqres_n <- t(array(c( 0, 0.1, barrierbase, 0.1*barrierbase, # down is starting resistance (-, c, n, cn), across is acquired pattern (-, c, n, cn) after novel regimen treatment
                                                                     0, 0, 0, 0.005, 
                                                                     0, 0, 0, 0.1, 
                                                                     0, 0, 0, 0), dim=c(4,4))); TRP$barrier$optimal$acqres_n[TRP$barrier$optimal$acqres_n>1] <- 1
    TRP$exclusions$minimal$eligibility <- 1- rep(0.11, 2)
    TRP$exclusions$optimal$eligibility <- 1- rep(0, 2)
    TRP$exclusions_hiv$minimal$eligibility[2] <- 0
    TRP$exclusions_hiv$optimal$eligibility[2] <- TRP$exclusions_hiv$optimal$eligibility[2] - 0.05
    
    TRP$tolerability$minimal$ltfurate_n <- values$ltfurate_sr
    TRP$tolerability$optimal$ltfurate_n <- values$ltfurate_sr - 0.03/values$months_n
  }
  
  if (target=="DR")
  {
    poor_n <- 0.18; TRP$effectiveness$minimal$poor_n <- c(poor_n + c(0, 0.13, 0.5), 1)
    poor_n <- 0.03; TRP$effectiveness$minimal$poor_n <- c(poor_n + c(0, 0.13, 0.5), 1)
    
    TRP$duration$minimal$months_n <- 18
    TRP$duration$optimal$months_n <- 6
    
    TRP$companion$minimal$cres <- rep(0.25, 2)
    TRP$companion$optimal$cres <- rep(0, 2)
    
    barrierbase <- 0.1; TRP$barrier$minimal$acqres_n <- t(array(c( 0, 0.1, barrierbase, 0.1*barrierbase, # down is starting resistance (-, c, n, cn), across is acquired pattern (-, c, n, cn) after novel regimen treatment
                                                                    0, 0, 0, 8*barrierbase, 
                                                                    0, 0, 0, 0.1, 
                                                                    0, 0, 0, 0), dim=c(4,4))); TRP$barrier$minimal$acqres_n[TRP$barrier$minimal$acqres_n>1] <- 1
    barrierbase <- 0.008; TRP$barrier$minimal$acqres_n <- t(array(c( 0, 0.1, barrierbase, 0.1*barrierbase, # down is starting resistance (-, c, n, cn), across is acquired pattern (-, c, n, cn) after novel regimen treatment
                                                                 0, 0, 0, 0.005, 
                                                                 0, 0, 0, 0.1, 
                                                                 0, 0, 0, 0), dim=c(4,4))); TRP$barrier$optimal$acqres_n[TRP$barrier$optimal$acqres_n>1] <- 1
    TRP$exclusions$minimal$eligibility <- 1- rep(0.11, 2)
    TRP$exclusions$optimal$eligibility <- 1- rep(0, 2)
    TRP$exclusions_hiv$minimal$eligibility[2] <- 0
    TRP$exclusions_hiv$optimal$eligibility[2] <- TRP$exclusions_hiv$optimal$eligibility[2] - 0.05
    
    TRP$tolerability$minimal$ltfurate_n <- values$ltfurate_sr
    TRP$tolerability$optimal$ltfurate_n <- values$ltfurate_sr - 0.06/values$months_n
  }
  
return(TRP)
}
  



set.values <- function()
{
  values <- list()
  values <- within(values, {
    selfcurerate <- 0.2; relapserate <- 1; latreduct <- 0.5 #self cure in hiv neg only
    hivrate <- 0.001; beta <- 8;  #don't matter here; will fit to target coprevalence
    mort <- c(0.012,0.033);  #by HIV status; 
    reactrate <- c(0.0015,0.03); rapidprog <- c(0.13,0.5); 
    transmissibility <- 0.7 #of DR strain(s); will LHS sample
    tbmort <- c(0.2,1) #by HIV status. LHSampled below for non-HIV. 
  
    months_s <- 6; months_r <- 18
    acqres_s <- 0.005; acqres_r <- 0
    poor_s <- c(0.06, 0.8); poor_r <- 0.24 #fraction with poor outcomes of either relapse or failure/tbdeath, of those who don't acquire resistance, for each relevant initial resistance pattern, (with the below relapsepoor fraction of the poor outcomes being relapses)
    relapsepoor <- 0.6
    dxrate <- c(0.67,0.9, # new  hiv- and hiv+, previously treat hiv- and hiv+
                      2,2.7)
    DSTrif <- c(0.1,1); # sampled later in LHS
    ltfurate_sr <- 0.01; # sampled later
    relapse246 <- c(7.5, 3, 1);   names(DSTrif) <- c("An","Ap"); names(relapse246) <- c("2mo","4mo","6mo") #extra relapse added by thirds of course completed (will interpolate between these) -- multiplicative by regimen efficacy
    
    availability <- 1 # Current version (editable within dxdt (nvary) scales up to this over 3 years 
    targetpop <- c(1,1); names(targetpop) <- c("DS", "DR") #will change for DR=(0,1), DS=(1,0), or panTB=(1,1) #actually set in samplenovel
    cres <- c(0.05,0.1); names(cres) <- c("rifs", "rifr") #baseline companion resistance prevalence among rif S and rif R; actually set in samplenovel.
    DSTnew <- c(0,1); names(DSTnew) <- c("c","n") #companion drugs, novel drug (doesn't depend on rif or retreatment status); actually set in samplenovel
    months_n <- 4 #actually set in samplenovel
    eligibility <- c(1,0.9)  #by HIV status #actually set in samplenovel
    ltfurate_n <- ltfurate_sr
    poor_n <- c(0.03 + c(0, 0.13, 0.5), 1) #for res to (0, c, n, cn), of those who don't acquire resistance #actually set in samplenovel
    acqres_n <- t(array(c( 0, 0.1, 0.05, 0.005, # down is starting resistance (-, c, n, cn), across is acquired pattern (-, c, n, cn) after novel regimen treatment #actually set in samplenovel
                           0, 0, 0.4, 0.04, 
                           0, 0, 0, 0.1, 
                           0, 0, 0, 0), dim=c(4,4)))
  })
  return(values)
}




samplepars <- function(whichparset)
{
  if (whichparset=="ds") samplepars <- t(array(c(
    "reactrate", 1, 0.0005,0.0025,
    "rapidprog", 1, 0.09,0.17, 
    "tbmort",1, 0.1,0.3,
    "poor_s",1, 0.01,0.11,
    "relapse246",2,  1.5, 4.5,
    "dxrate",1,  0.3,1,
    "ltfurate_sr", 1, 0.005, 0.015
  ), dim=c(4,7) )) else
    
    if (whichparset=="dr") samplepars <- t(array(c(
      "acqres_s", 1, 0.001, 0.009,
      "poor_r", 1, 0.13, 0.35, #  this is of those who are susceptible, remain on treatment, and don't acquire additional resistance. i've estimated 25% (see appendix for sources) and am using WHO global report for range. Globally, about 2/3 of those with known outcomes are successes rather than failure or death. But some of that may be due to XDR, etc. Minimum poor assuming 20% resistance-related failures and no relapses is ~13%, and maximum assuming no resistance-related failures and assuming relapses=failures is 35%. 
      "transmissibility", 1, 0.5, 0.9,
      "DSTrif", 1, 0.1, 0.5, # make DST time varying over five years, with these as max of 5-year linear increase. Then hold constant once novel regimen introduced.
      "DSTrif", 2, 0.6, 1 # 
    ), dim=c(4, 5) )) else 
      
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



screendrout <- function(drout_filename=paste0("DRcalibration_",tag,".csv"))
{
  drout <- read.csv(file = drout_filename, header = TRUE) #saved results from dr sampling runs to time 0
  screened <- drout[ drout[,"rrinc"]/drout[,"inc"] > 1/3*drout[,"targetdr"] & drout[,"rrinc"]/drout[,"inc"] < 3*drout[,"targetdr"] ]  #within 3fold if rr incident fraction target
  return(screened)
}


evaltrp <- function(genericvalues, drsetup, drout, ids, idr, targetpt="DS", DST=FALSE)
{
  
  if(missing(genericvalues)) {genericvalues <- readRDS(paste0("genericvalues_",tag,".RDS"))} # source of parameters that have fixed values throughout (as saved at start of sampling)
  if(missing(drsetup)) {drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)}
  if(missing(drout)) {drout <- read.csv(file = paste0("DRcalibration_", tag, ".csv"), header = TRUE)} #includes i,i,targetepi,beta, hivrate, sampledpars(ds/dr), finalstate
  
  rows<-1;  if (missing(ids) || missing(idr)) {rows <- 1:nrow(drout)} else {rows <- 1:nrow(drout)[(drout[,"ids"] %in% ids) && (drout[,"idr"] %in% idr)]}
  
  novelsetup <- setup.model(DRera = TRUE, treatSL = TRUE, treatnovel = TRUE)
  
  allsampledpars <- c("beta", "hivrate", unique(rbind(samplepars("ds"), samplepars("dr"))[,1]))
  nsampledpars <- 0; for (n in 1:length(allsampledpars)) { nsampledpars <- nsampledpars + length(genericvalues[[allsampledpars[n]]]) }
  
  write(c("inew", "ids","idr","targetepi","targetpt","DST", "varied element", "level", tallynames),  #create header for out file (make sure this is updates for any edits to iter or to later write function)
        file=paste0("TRPoutput_", targetpt,DST,"_",tag,".csv"), sep=",", ncol=8+length(tallynames))
  
  for (inew in rows)
  {
    iter <- unlist(c(unlist(drout[inew,c("ids", "idr", "targetprev","targetcoprev","targetdr")]), targetpt, DST)) #will include these labels as part of returned output
    
    valuevect <- drout[inew, 5+(1:nsampledpars)] 
    v <- 0 #initialize at start of valuevect
    for (pname in allsampledpars)
    {
      genericvalues[[pname]] <- unlist(valuevect[v+(1:length(genericvalues[[pname]]))])
      v <- v + length(genericvalues[[pname]]) # move forward to start of next par vector in sampled values
    }    
    
    state <- drout[inew, 5+nsampledpars + (1:length(drsetup$statenames))] #5=isimds, isimdr, targetepi, beta, and hivrate
    
    # now we've pulled everything back out from drout. next, we need to set up for including novel regimen, e.g. expand the novel to include novel treatment and resistance. 
    newstate <- numeric(length(novelsetup$statenames)); names(newstate) <- novelsetup$statenames
    newstate[drsetup$statenames] <- unlist(state)
    
    # and now sample the TRP -- generating a baseline result and a minimal and optimal for each TRP element
    TRP <- samplenovel(genericvalues, targetpt, DST) #a list (by element) of lists (by level) of values
    
    for (vary in names(TRP))
    {  for (level in names(TRP[[vary]]))
      {
        valueset <- TRP[[vary]][[level]]
        
        # revise state to assign companion resistance to specified fractions of 0 and r's
        s_cr <- valueset$cres["rifs"]; r_cr <- valueset$cres["rifr"]
        novelstate <- newstate * rep( rep(c(1-s_cr, s_cr, 1-s_cr, s_cr, 1-r_cr, r_cr, 1-r_cr, r_cr), each=length(novelsetup$statenames)/16), 2)
        
        parset <- create.pars(setup = novelsetup, values = valueset, T, T, T)
        
        ## implement novel regimen with the given TRP, and record yearly state and stats for ten years
        outset <- ode(y=unlist(novelstate), times=seq(0, 10, by=0.1), func=dxdt, parms=parset$fullpars, do.tally=TRUE, method=lsodes)
        
        TRP[[vary]][[level]]$output <- outset
        
        write(c( iter, vary, level, outset[11,1+length(novelstate)+(1:length(tallynames))]), 
              file=paste0("TRPoutput_", targetpt,DST,"_",tag,".csv"), sep=",", append=TRUE)
      
      }
    } 
    
  }
  return(TRP)
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
      
      dxrate <-  t(array(dxrate, dim=c(2,2))); colnames(dxrate) <- c("An","Ap"); rownames(dxrate) <- Hnames
      
      names(reactrate) <- names(rapidprog) <- names(mort) <- names(tbmort) <- Hnames  
      
      transmissibility <- c(1,rep(transmissibility, length(Rnames)-1)); names(transmissibility) <- Rnames #need to decide details of fitness costs and sampling
      
      # then adjust for whether regimens are in use (i.e.if not using FL, no diagnosis and treatment; if not using SL, no rif diagnosis (and therefore no second-line treatment); 
      # and if not using novel regimen (which can be used in absence of novel DST in the model), availability of novel regimen is 0.
      if (DRera==FALSE) acqres_s <- 0; if (treatSL==FALSE) DSTrif <- c(0,0); if (treatnovel==FALSE) availability <- 0 
      
      names(eligibility) <- Hnames
      ltfurate <- c(ltfurate_sr, ltfurate_sr, ltfurate_n); names(ltfurate) <- regimens
      
      
      #translate above parameters into outcome arrays (shouldn't need to edit this part unless model structure changes):
      poormat <- array(rbind(
        c(rep(poor_s[1], max(1,4*treatnovel)), rep(poor_s[2], DRera*max(1,4*treatnovel))),
        poor_r,
        poor_n), dim=c(3,length(Rnames))); dimnames(poormat) <- list(regimens, Rnames) 
      
      failmat <- (1-relapsepoor)*poormat
      relapsemat <- relapsepoor*poormat
      
      acqresmat <- array(0,dim=c(length(Rnames),length(Rnames),length(regimens))); dimnames(acqresmat)=list(Rnames, Rnames, regimens) # from old resistance (down) to new resistance (across), by regimen 
      acqresmat[grep("^R[0cn]",Rnames), grep("^Rr",Rnames),"s"] <-  acqres_s * diag(max(1, 4*treatnovel))
      acqresmat[,,"r"] <- array(0,dim=c(length(Rnames),length(Rnames))) # currently no acqres with regimen r
      if (treatnovel) { acqresmat[,,"n"] <- rbind(cbind(acqres_n,0,0,0,0), cbind(0,0,0,0, acqres_n))[1:length(Rnames), 1:length(Rnames)] }
      
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
     
    # self cure (HIV neg only, and all to R0 resistance), TB reactivation, mortality (background and TB-related, placed as ins on diagonal for now), and relapse:
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
      }
    }
    

    ## TB diagnosis and treatment initiation happen in dxdt (to allow time-varying treatment availability)
    
    # once on effective treatment, progress through treatment to either relapse or cure (or to active if lost to follow up during first 2 months), including a monthly rate of loss to follow up 
    # this function determines whether relapse or cure depending on effectiveness of regimen and fraction completed. Assuming losses occur in the middle of each time period on average.
    relapsefracs <- function(period) # use relapse %s at 2, 4, and 6 months (defined in pars) for a 6 month regimen, and interpolate linearly 2--4 and 4--6 to get relapse % after any fraction of treatment course
    {
      fraction_completed <- (sum(Tperiod.lengths[1:period-1])+1/2*Tperiod.lengths[period])/durations
      relapse <- array(0, dim=c(length(Rnames),3)); dimnames(relapse) <- list(Rnames, regimens)
      relapse[,fraction_completed < 1/3] <- 1 
      relapse[,fraction_completed >= 1/3 & fraction_completed < 2/3] <- 
        (t((relapsemat) * (relapse246[2] + (relapse246[1]-relapse246[2])*(2/3-fraction_completed))))[, fraction_completed >= 1/3 & fraction_completed < 2/3]
      relapse[,fraction_completed >= 2/3 & fraction_completed < 1] <- 
        (t((relapsemat) * (relapse246[3] + (relapse246[3]-relapse246[3])*(1-fraction_completed))))[, fraction_completed >= 2/3 & fraction_completed < 1]
      relapse[,fraction_completed >= 1] <- t(relapsemat)[, fraction_completed >= 1]
      relapse[relapse>1] <- 1 #if sum exceeds 100%, set to 100% relapse
      return(relapse) 
    }
    

    for (jh in 1:length(Hnames))
    {
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
# for rvary and nvary, t needs to be defined with respect to 0 as the end of rifDST scaleup and the time of novel regimen implementation
# output is ode output, along with outcomes (as defined within function) at each time point
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
    if (length(Rnames)==1) { FOI = sum(statemat[c("An","Ap","Ti"),,]) * transmissibility * beta / sum(statemat) 
      } else  { FOI <- apply(statemat[c("An","Ap","Ti"),,], 2, sum) * transmissibility * beta / sum(statemat) }# FOI by strain
    
    # infection
      if (length(Rnames)>1 && sum(statemat[c("S","C"), -1, ]) > 0.0001 ) { stop("Error: Some susceptibles have drug resistance and won't be included in infection events.") }
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
          #acquired resistance moves to pending relapse for *new* strain
          mat[it, "R", , , jh, jh]  <- 
            mat[it, "R", , , jh, jh] + dxrate[jh, it] * startmat[it,Rnames,nreg] * acqresmat[,,nreg]
          
          #of those who don't acquire resistance (1-rowSums(acqresmat)), will split between Ti (failmat) and T1 for the selected (startmat) regimen
          if (length(Rnames)>1)
          {
            diag(mat[it, "Ti", , , jh, jh]) <- 
              diag(mat[it, "Ti", , , jh, jh]) + dxrate[jh, it] * startmat[it, Rnames ,nreg] * (1-rowSums(acqresmat[,,nreg]))*failmat[nreg, ] #failures go to ineffective treatment (with ongoing infectiousness and increased mortality risk)
            diag(mat[it, grep("T[srn]1",Tnames)[nreg], , , jh, jh]) <- 
              diag(mat[it, grep("T[srn]1",Tnames)[nreg], , , jh, jh]) + dxrate[jh, it] * startmat[it, Rnames, nreg] * (1-rowSums(acqresmat[,,nreg]))*(1 - failmat[nreg, ]) #for the remainder, treatment is initially effective (i.e. outcome will be cure vs relapse if they complete at least 2 months)
          } 
          else #if only the standard regimen is an option, diag and rowSums functions cause errors
          {
            mat[it, "Ti", 1,1, jh, jh] <- 
              mat[it, "Ti", 1,1, jh, jh] + dxrate[jh, it] * startmat[it, Rnames ,nreg] * (1-(acqresmat[,,nreg]))*failmat[nreg,Rnames ] #failures go to ineffective treatment (with ongoing infectiousness and increased mortality risk)
            mat[it, grep("T[srn]1",Tnames)[nreg], , , jh, jh] <- 
              mat[it, grep("T[srn]1",Tnames)[nreg], , , jh, jh] +  dxrate[jh, it] * startmat[it, Rnames, nreg] * (1-(acqresmat[,,nreg]))*(1 - failmat[nreg, Rnames]) #for the remainder, treatment is initially effective (i.e. outcome will be cure vs relapse if they complete at least 2 months)            
          }
        }
      }
    }
    
  
    ########### transform mat
    newmat <- aperm(mat, c(1,3,5,2,4,6)) # this will translate into a 2d matrix more easily
    squaremat <- array(newmat, dim=c(length(Tnames)*length(Rnames)*length(Hnames), length(Tnames)*length(Rnames)*length(Hnames))) #2d state1 to state2 transmition matrix
    dimnames(squaremat) <- list(statenames, statenames)
    
    
    ########## tally how much each current state contributes to outcomes of interest (then will multiply by state and append to output)
    
    outcomes <- c("prev", "inc", "rrinc", "rronsets", "panronsets", "relapses", "tbdeaths", "rrdeaths", "dxs", "rDSTs", "nDSTs", "rxtime_s", "rxtime_r", "rxtime_n")
    tally <- array(0,dim=c(length(Tnames)*length(Rnames)*length(Hnames), length(outcomes))); dimnames(tally) <- list(statenames, outcomes)
    
    if (do.tally==TRUE)
    {
      tally[c(grep("^A", statenames), grep("^T", statenames)), "prev"] <- rep(1, length(c(grep("^A", statenames), grep("^T", statenames))))
      
      tally[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames)),"inc"] <- #doens't include relapses (mostly early), but does include reinfections (mostly late) - i.3. this is a count of new infecitons, whether or not they are recognized as such
        apply(squaremat[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames)),  grep("^A",statenames)], 1, sum)
    
      tally[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames)),"rrinc"] <- #also only counts new infections with rr strain
        apply(squaremat[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames)),  grep("^A.[.]Rr",statenames)], 1, sum)
    
      tally[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames), grep("^R",statenames)),"rronsets"] <- #this alternative also includes resistance relapses (including resistance acquisitions with relapse)
        apply(squaremat[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames),grep("^R",statenames)),  grep("^A.[.]Rr",statenames)], 1, sum)
      
      tally[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames), grep("^R",statenames)),"panronsets"] <-
        apply(squaremat[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames),grep("^R",statenames)),  grep("^A.[.]R(r|rc|rcn|rn|cn|n)",statenames)], 1, sum)

      tally[grep("^R",statenames),"relapses"] <- 
        apply(squaremat[grep("^R",statenames),  grep("^A",statenames)], 1, sum)
      
      tally[c(grep("^A", statenames), grep("^Ti", statenames)), "tbdeaths"] <- rep(tbmort, each=(length(c(grep("^A", statenames), grep("^Ti", statenames)))/2))
      
      tally[c(grep("^A.[.]Rr", statenames), grep("^Ti[.]Rr", statenames)), "rrdeaths"] <- rep(tbmort, each=(length(c(grep("^A.[.]Rr", statenames), grep("^Ti[.]Rr", statenames)))/2))
      
      tally[grep("^A", statenames),"dxs"] <- apply( squaremat[ grep("^A", statenames), c(grep("^T", statenames), grep("^R", statenames)) ], 1, sum)
      
      tally[grep("^A", statenames),"rDSTs"] <- 
        apply( squaremat[ grep("^A", statenames), c(grep("^T",statenames),grep("^R",statenames)) ] * DSTrif_t, 1, sum) #alternates between An and Ap, so alternate DST coverage
               
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
  }
 )
}



#advance forward a bit:
## stateplus includes outcomes tracking variables
advance <- function(state, t0, addedt=1, rvary, nvary, reportsteps=1, fullpars) 
{
  o <- ode(state, seq(t0,t0+addedt,by=0.1), dxdt, parms=fullpars, rvary=rvary, nvary=nvary, do.tally=TRUE, method=lsodes)
  return(o[(10*(addedt)+1 - reportsteps):(10*(addedt)+1),])
}

  
# run to equilibrium without drugs; will need to adjust Rnames and resistance-related parameters, and remake mat
equilib <- function(state, pars, tol=0.2)
{
  if(missing(pars)) pars <- create.pars(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
  
  if(missing(state)) { state <- with(pars$fullpars,  c(90000, 9900, 100, rep(0,length(statenames)-3))); names(state) <- pars$fullpars$statenames }
  
  statex <- ode(state, seq(0,20,by=0.1), func=dxdt, pars$fullpars, do.tally=TRUE, method=lsodes)[200:201,]; totaltime <- 20
  
  log <- c(0, state, rep(0, ncol(statex)-1-length(state))); names(log) <- colnames(statex)
  log <- rbind(log, statex[2,])
  
  while (max (abs(statex[2,-1]-statex[1,-1]))>tol)
  {
    statex <- advance(state=statex[2,2:(length(state)+1)], t0=totaltime, addedt=10, rvary=FALSE, nvary=FALSE, fullpars=pars$fullpars)
    totaltime <- totaltime + 10
    log <- rbind(log, statex[2,])
  }
  return(list("log"=log, "pars"=pars))
}
  
# plotlog <- function(log)
# {
#   par(mfrow=c(2,2), mar=c(2,4,1,1))
#   plot(log[,1], log[,"inc"], type='l') 
#   plot(log[,1], log[,"prev"]) 
#   plot(log[,1], log[,"tbdeaths"]) 
#   plot(log[,1], log[,"dxs"]) 
#   
# }

# adddr<- function(eqb, acqres=0.005, DSTrif = c(0,0), rtrans=0.7, addedt=10)
# {
#   if(missing(eqb)) eqb <- equilib()
#   
#   oldstatenames <- eqb$pars$fullpars$statenames
#   oldstate <- with(eqb,log[nrow(log),2:(length(oldstatenames)+1)])
#   
#   oldvalues <- eqb$pars$values
#   newvalues <- oldvalues; newvalues$acqres_s <- acqres; newvalues$DSTrif <- DSTrif; newvalues$transmissibility <- rtrans
#   
#   if (sum(DSTrif)>0) treatSL <- TRUE else treatSL <- FALSE
#   newpars <- create.pars(values=newvalues, DRera=TRUE, treatSL=treatSL, treatnovel=FALSE)
#   
#   newstate <- numeric(length(newpars$fullpars$statenames)); names(newstate) <- newpars$fullpars$statenames
#   newstate[oldstatenames] <- oldstate
#   
#   drode <- advance(newstate, t0=0, dxdt, addedt, reportsteps = 10*addedt, fullpars=newpars$fullpars) 
#   
#   par(mfrow=c(2,2), mar=c(2,4,1,1))
#   plot(drode[,"time"], drode[,"inc"])
#   plot(drode[,"time"], drode[,"rrinc"])
#   plot(drode[,"time"], drode[,"rrdeaths"])
#   plot(drode[,"time"], drode[,"rDSTs"])
#   
#   return(list("drode"=drode, "pars"=newpars))
# }

