library("deSolve")
library("stats") #contains optim()
library("lhs")
library("akima")
library("stringr")

targetepis <- list("Brazil"=c(52, 0.17, 0.014), "India"=c(195, 0.04, 0.022), "Philippines"=c(417, 0.002, 0.02), "SouthAfrica"=c(696, 0.61, 0.018)) # tb prev, HIV coprev, rrinc/inc

# will provide a list of "values" inputs for the function that implements novel regimen: 
# each labeled by the targetpop, the DST use, the TRP element varied (with "none" as one option), and whether the varied TRP element is minimal, intermediate, or optimal

set.novelvalues <- function()
{
  selections <- list()
  selections$poor_n <- array(c(0.06, 0.03, 0, 
                               0.24, 0.06, 0.03), dim=c(3,2))
  selections$months_n <- array(c(6, 4, 2, 
                                 20, 9, 6), dim=c(3,2))
  selections$cres <- array(c(0.1,0.1,  0.03,0.03,  0,0, 
                             0.15,0.15,  0.05,0.05,  0,0), dim=c(2,3,2))
  selections$barrierbase <- array(c(0.05, 0.008, 0,
                                    0.1, 0.05, 0.008), dim=c(3,2))
  selections$eligibility <- array(c(0.89,0,  0.95,0.9,  1,1,
                                    0.89,0,  0.95,0.9,  1,1), dim=c(2,3,2))
  # !! can add specifics e.g. hiv or pediatric low/med/high as separate selection item
  selections$ltfu_reduction <- array(c(0, 0.015, 0.03,
                                       0, 0.03, 0.06), dim=c(3,2))
  elementnames <- c("efficacy", "duration", "companion", "barrier", "exclusions", "tolerability")
  return(list("selections"=selections, "elementnames"=elementnames))
}


set.values <- function()
{
  values <- list()
  values <- within(values, {
    cal <- list(); cal <- within(cal, {
      hivrate <- 0.01; beta <- 10;  #will fit to target coprevalence
    })
    varied_ds <- list(); varied_ds <- within(varied_ds, {
      selfcurerate <- 0.2; relapserate <- 1; latreduct <- 0.6 #self cure in hiv neg only
      mort <- c(0.012,0.033);  #by HIV status; 
      reactrate <- c(0.0015,0.03); rapidprog <- c(0.13,0.5); 
      tbmort <- c(0.2,1) #by HIV status. 
      poor_s_rifs <- 0.06
      relapsepoor <- 0.67
      dxrate <- c(0.6,0.9, # new and rerx  hiv-, and new and rerx hiv+
                  2, 3)
      ltfurate_sr <- 0.01 #monthly rate
      initialloss_s <- 0.15
      relapse24 <- c(7.5, 3)
      transcost_n <- 0.3
    })
    varied_dr <- list(); varied_dr <- within(varied_dr, { #includes those vary for novel regimen outside of trp
      nonpoor_s_rifr <- 0.2
      acqres_s <- 0.008; 
      poor_r <- 0.24 #fraction with poor outcomes of either relapse or failure/tbdeath, of those who don't acquire resistance, for each relevant initial resistance pattern, (with the below relapsepoor fraction of the poor outcomes being relapses)
      transcost_rif <- 0.3 #redeuction (from 1) in fitness for RifR strain(s)
      DSTrif_n <- 0.3
      noDSTrif_p <- 0.3
      poor_n_cnr <- c(0.14,0.2)
      acqres_candn <- 0.1# probably of acquiring c resistance if already n resistant (in same or subsequent treatment course)
      acqres_nifc <- 6 # increase in probably of acquiring n resistance if already c resistant (in same treatment course)
      
    })
    fixed <- list(); fixed <- within(fixed, { #includes those that will vary with novel regimen within trp
      months_s <- 6; months_r <- 20
      acqres_r <- 0
      availability <- 0.8 # Current version (editable within dxdt (nvary) scales up to this over 3 years 
      targetpop <- c(1,0) #will change for DR=(0,1), DS=(1,0), or panTB=(1,1) #actually set in samplenovel
      cres <- c(0.05,0.1) #baseline companion resistance prevalence among rif S and rif R; actually set in samplenovel.
      DSTnew <- c(1,1) #companion drugs, novel drug (doesn't depend on rif or retreatment status); actually set in samplenovel
      months_n <- 4 #actually set in samplenovel
      eligibility <- c(1,1)  #by HIV status #actually set in samplenovel
      ltfurate_n <- varied_ds$ltfurate_sr
      poor_n <- 0.03 #will set value in samplenovel
      acqres_n <- rep(0, times=16) #will set values in samplenovel, based in part on varied_dr values above
    })
  })
  return(values)
}


#returns a set of values for a specified starting values and trp element combination
sampleTRP <- function(mergedvalues, targetpt="DS", DST="DSTall", optimals=NA, minimals=NA, HIV="both")
{
  selections <- set.novelvalues()$selections
  elementnames <- set.novelvalues()$elementnames
  
  if (DST=="DSTall") {mergedvalues$DSTnew[1:2] <- c(1,1)} else {mergedvalues$DSTnew[1:2] <- c(0,0)}
  
  if (targetpt == "DS") {mergedvalues$targetpop <- c(1,0); whichcol<-1} else {mergedvalues$targetpop <- c(0,1); whichcol<-2 }
  whichrow <- rep(2, length(elementnames)); whichrow[which(elementnames %in% optimals)] <- 3; whichrow[which(elementnames %in% minimals)] <- 1
  
  mergedvalues$poor_n <- selections$poor_n[whichrow[which(elementnames=="efficacy")], whichcol]
  mergedvalues$months_n <- selections$months_n[whichrow[which(elementnames=="duration")], whichcol]
  mergedvalues$cres[1:2] <- selections$cres[ , whichrow[which(elementnames=="companion")], whichcol]
  barrierbase <- selections$barrierbase[whichrow[which(elementnames=="barrier")], whichcol];
  mergedvalues$acqres_n <- t(array(c( 0, 0, (1-mergedvalues$acqres_candn)*barrierbase, mergedvalues$acqres_candn*barrierbase, # down is starting resistance (-, c, n, cn), across is acquired pattern (-, c, n, cn) after novel regimen treatment
                                      0, 0, 0, mergedvalues$acqres_nifc*barrierbase, 
                                      0, 0, 0, mergedvalues$acqres_candn, 
                                      0, 0, 0, 0), dim=c(4,4))); if (barrierbase==0) mergedvalues$acqres_n[2,4] <- 0.005; mergedvalues$acqres_n[mergedvalues$acqres_n>1] <- 1
  mergedvalues$eligibility[1:2] <- selections$eligibility[ , whichrow[which(elementnames=="exclusions")], whichcol]
  if (HIV=="HIV") mergedvalues$eligibility[1] <- 1
  if (HIV=="nonHIV") mergedvalues$eligibility[2] <- mergedvalues$eligibility[1]
  ltfu_reduction <- selections$ltfu_reduction[whichrow[which(elementnames=="tolerability")], whichcol]; 
  mergedvalues$ltfurate_n <- max(0, mergedvalues$ltfurate_sr - ltfu_reduction/mergedvalues$months_n)
  
  return(mergedvalues) # a set of edited single-list values for use in create.pars
}

# evaluates 10-year impact for optimal and minimal for each TRP element, with others at intermediate ("typical") level, and full trajectory for variation in all TRP parameters together
evaltrp <- function(genericvalues, drsetup, drout, ids, idr, rows, targetpt="DS", DST="DSTall", tag=currenttag,location="",rDSTall=FALSE) # uses merged but not unlisted values
{
  if(missing(genericvalues)) {genericvalues <- readRDS(paste0("genericvalues_",tag,".RDS"))} # source of parameters that will have constant values (as saved at start of sampling)
  if(length(genericvalues[[1]][[1]]>1)) { genericvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]])) } #merge into single list if needed
  if(missing(drsetup)) {drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)}
  if(missing(drout)) {drout <- read.csv(file = paste0("DRcalibration_", tag, ".csv"), header = TRUE)} #includes i,i,targetepi,beta, hivrate, sampledpars(ds/dr), finalstate
  
  if(missing(rows))  if (missing(ids) || missing(idr)) {rows <- 1:nrow(drout)} else {rows <- (1:nrow(drout))[(drout[,"ids"] %in% ids) & (drout[,"idr"] %in% idr)]}
  novelsetup <- setup.model(DRera = TRUE, treatSL = TRUE, treatnovel = TRUE)
  elementnames <- c("all", set.novelvalues()$elementnames)
  
  stateheader <- c("inew", "ids","idr","targetprev","targetcoprev", "targetdr", "targetpt","DST", "vary","level","time", novelsetup$statenames)
  
  wideheader <- c("inew", "ids","idr","targetprev","targetcoprev", "targetdr", "targetpt","DST", names(unlist(genericvalues)))
  wideheader <- append(wideheader, paste0(rep(tallynames,times=11*3),rep(rep(0:10, each=length(tallynames)), times=3), 
                                          rep(c("allminimal", "allintermediate","alloptimal"), each=11*length(tallynames))))
  for (i in 2:length(elementnames)) wideheader <- append(wideheader, 
                                                         paste0( rep(tallynames, times=2*11), 
                                                                 rep( rep(0:10, each=length(tallynames)), 2),
                                                                 rep(elementnames[i], each=22*length(tallynames) ),
                                                                 rep( c("minimal","optimal"), each=11*length(tallynames) ) ) )
  
  
  if(!file.exists(paste0(location,"TRPstateoutput_", targetpt,DST,"_",tag,".csv"))) { write(stateheader,  file=paste0(location,"TRPstateoutput_", targetpt,DST,"_",tag,".csv"), sep=",", ncol=length(stateheader)) }
  if(!file.exists(paste0(location,"TRPwideoutput_", targetpt,DST,"_",tag,".csv"))) { write(wideheader,  file=paste0(location,"TRPwideoutput_", targetpt,DST,"_",tag,".csv"), sep=",", ncol=length(wideheader)) }
  
  for (inew in rows)
  {
    iter <- unlist(c(inew,unlist(drout[inew,c("ids", "idr", "targetprev","targetcoprev","targetdr")]), targetpt, DST)) #will include these labels as part of returned output
    
    
    valuevect <- unlist(drout[inew, 5+(1:length(unlist(genericvalues)))]) ; names(valuevect) <- names(unlist(genericvalues))
    if (rDSTall==TRUE) { valuevect[which(names(valuevect)=="DSTrif_n.varied_dr.DSTrif_n")] <- 1; valuevect[which(names(valuevect)=="noDSTrif_p.varied_dr.noDSTrif_p")] <- 0 }
    v <- 0; for (pname in names(genericvalues))
    { genericvalues[[pname]] <- valuevect[v+(1:length(genericvalues[[pname]]))]; #names(genericvalues[[pname]]) <- names()
      v <- v + length(genericvalues[[pname]]) # move forward to start of next par vector in sampled values
    }    
    
    state <- drout[inew,drsetup$statenames]
    newstate <- numeric(length(novelsetup$statenames)); names(newstate) <- novelsetup$statenames
    newstate[drsetup$statenames] <- unlist(state)
    
    iresult <- numeric(0) #will become full row of results for this dr row (inew), for all elements and levels
    for (vary in elementnames)
    {
      print(paste0("Evaluating TRP variation in ",vary,  ",  Simulation #", inew," of ",nrow(drout)," for ",targetpt,DST,currenttag))
      valueset <- list()
      if (vary=="all")
      { valueset$m <- sampleTRP(mergedvalues = genericvalues, targetpt = targetpt, DST = DST, minimals=elementnames[-1])
        valueset$i <- sampleTRP(mergedvalues = genericvalues, targetpt = targetpt, DST = DST)
        valueset$o <- sampleTRP(mergedvalues = genericvalues, targetpt = targetpt, DST = DST, optimals=elementnames[-1])
      } else 
      { valueset$m <- sampleTRP(mergedvalues = genericvalues, targetpt = targetpt, DST = DST, minimals=vary)
        valueset$o <- sampleTRP(mergedvalues = genericvalues, targetpt = targetpt, DST = DST, optimals=vary)
      }
      for (level in names(valueset))
      {
        s_cr <- valueset[[level]]$cres[1]; r_cr <- valueset[[level]]$cres[2]
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
        
        parset <- create.pars(setup = novelsetup, values = valueset[[level]], T, T, T)
        
        outset <- ode(y=unlist(novelstate), times=0:10, func=dxdt, parms=parset$fullpars, do.tally=TRUE, method="adams")
        
        write(t(cbind(t(array(rep(c(iter,vary,level), times=11),dim=c(length(iter)+2,11))), outset[,c("time",novelsetup$statenames)])), 
                                file=paste0(location,"TRPstateoutput_",targetpt,DST,"_",tag,".csv"), sep=",", append=TRUE, ncol=length(stateheader))
          
        iresult <- append(iresult, as.vector(t(outset[,tallynames]))) 
          
      }
    } 
    
    write(unlist(c(iter, valuevect, iresult)), file=paste0(location,"TRPwideoutput_", targetpt,DST,"_",tag,".csv"), sep=",", append=TRUE, ncol=length(wideheader))
    
  }
  return(0)
}



sample.values <- function(values, whichparset, LHS, isim)
{
  midvalues <- set.values()
  j <- 0
  for (n in 1:length(values[[whichparset]]))
  {
    values[[whichparset]][[names(values[[whichparset]])[n]]] <- 
      midvalues[[whichparset]][[names(values[[whichparset]])[n]]] * (0.5  + LHS[isim, j + (1:length(unlist(values[[whichparset]][n])))]) # = +/- up to 50% the midvalue of each in whichparset
    j <- j + length(unlist(values[[whichparset]][n]))
  }
  
  return(values)
}



screendrout <- function(drout_filename=paste0("../scratch/DRcalibration_",currenttag,".csv"), tolerance=3, usethistargetdr)
{
  drout <- read.csv(file = drout_filename, header = TRUE) #saved results from dr sampling runs to time 0
  if(missing(usethistargetdr))
  {
    screened <- drout[ drout[,"rrinc"]/drout[,"inc"] > 1/tolerance*drout[,"targetdr"] & drout[,"rrinc"]/drout[,"inc"] < tolerance*drout[,"targetdr"], ]  #within 3fold if rr incident fraction target
  } else screened <- drout[ drout[,"rrinc"]/drout[,"inc"] > 1/tolerance*usethistargetdr & drout[,"rrinc"]/drout[,"inc"] < tolerance*usethistargetdr, ]  #within 3fold if rr incident fraction target
  return(screened)
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
  N <- 8; lengths <- c(1,1,1,1,2,3,3,8)/12
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
  if (length(values)==4) pars <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]])) else pars <- values
  
  with(setup, {
    pars <- within(pars, {
      
      dxrate <-  array(unlist(dxrate), dim=c(2,2)); rownames(dxrate) <- c("An","Ap"); colnames(dxrate) <- Hnames
      
      names(reactrate) <- names(rapidprog) <- names(mort) <- names(tbmort) <- Hnames  
      
      transmissibility <- rep(1,length(Rnames)); 
      transmissibility[grep("r", Rnames)] <- transmissibility[grep("r", Rnames)]*(1-transcost_rif);
      transmissibility[grep("n", Rnames)] <- transmissibility[grep("n", Rnames)]*(1-transcost_n);
      
      relapse246 <- append(relapse24, 1)
      poor_s <- c(poor_s_rifs, 1-nonpoor_s_rifr)
      DSTrif <- c(DSTrif_n, 1-noDSTrif_p)
      initialloss <- rep(initialloss_s, 3)
      
      # then adjust for whether regimens are in use (i.e.if not using FL, no diagnosis and treatment; if not using SL, no rif diagnosis (and therefore no second-line treatment); 
      # and if not using novel regimen (which can be used in absence of novel DST in the model), availability of novel regimen is 0.
      if (DRera==FALSE) acqres_s <- 0; if (treatSL==FALSE) DSTrif <- c(0,0); if (treatnovel==FALSE) availability <- 0 
      
      names(eligibility) <- Hnames
      ltfurate <- c(ltfurate_sr, ltfurate_sr, ltfurate_n); names(ltfurate) <- regimens
      
      
      #translate above parameters into outcome arrays (shouldn't need to edit this part unless model structure changes):
      poormat <- array(rbind(
        c(rep(poor_s_rifs, max(1,4*treatnovel)), rep((1-nonpoor_s_rifr), DRera*max(1,4*treatnovel))),
        poor_r,
        c(poor_n + c(0, poor_n_cnr), 1)), dim=c(3,length(Rnames))); dimnames(poormat) <- list(regimens, Rnames) 
      
      failmat <- (1-relapsepoor)*poormat
      relapsemat <- relapsepoor*poormat
      
      acqresmat <- array(0,dim=c(length(Rnames),length(Rnames),length(regimens))); dimnames(acqresmat)=list(Rnames, Rnames, regimens) # from old resistance (down) to new resistance (across), by regimen 
      acqresmat[grep("^R[0cn]",Rnames), grep("^Rr",Rnames),"s"] <-  acqres_s * diag(max(1, 4*treatnovel))
      acqresmat[,,"r"] <- array(0,dim=c(length(Rnames),length(Rnames))) # currently no acqres with regimen r
      if (treatnovel==TRUE) {acqresmat[1:4,1:4,"n"] <- 
                               acqresmat[5:8,5:8,"n"] <- array(unlist(acqres_n), dim=c(4,4)) }
      acqresmat[acqresmat>1] <- 1
      
      
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
    
    # once on effective treatment, progress through treatment to either relapse or cure, including a monthly rate of loss to follow up 
    # this function determines whether relapse or cure depending on efficacy of regimen and fraction completed. Assuming losses occur in the middle of each time period on average.
    relapsefracs <- function(period) # use relapse %s at 2, 4, and 6 months (defined in pars) for a 6 month regimen, and interpolate linearly 2--4 and 4--6 to get relapse % after any fraction of treatment course
    {
      fraction_completed <- (sum(Tperiod.lengths[(1:length(Tperiod.lengths))<period])+1/2*Tperiod.lengths[period])/durations
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
# fullpars argument is the fullpars component of create.pars (i.e. list of params incl mat, but doesn't include list(s) of values)
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
    
    outcomes <- c("prev","rprev","cprev","nprev","cnprev","novelprev","rnovelprev", "inc", "rrinc", "rronsets", "panronsets", "relapses", "tbdeaths", "rrdeaths", "hivtoo", "dxs", "rDSTs", "nDSTs", "rxtime_s", "rxtime_r", "rxtime_n")
    tally <- array(0,dim=c(length(Tnames)*length(Rnames)*length(Hnames), length(outcomes))); dimnames(tally) <- list(statenames, outcomes)
    
    if (do.tally==TRUE)
    {
      tally[c(grep("^A", statenames), grep("^T", statenames)), "prev"] <- 1
      
      tally[c(grep("^A.[.]Rr", statenames), grep("^Ti[.]Rr", statenames)), "rprev"] <- 1
      
      tally[c(grep("^A.[.]R(rc|c)", statenames), grep("^Ti[.]R(rc|c)", statenames)), "cprev"] <- 1
      
      tally[c(grep("^A.[.]R(rcn|cn|n)", statenames), grep("^Ti[.]R(rcn|cn|n)", statenames)), "nprev"] <- 1
      
      tally[c(grep("^A.[.]R(rcn|cn)", statenames), grep("^Ti[.]R(rcn|cn)", statenames)), "cnprev"] <- 1
      
      tally[c(grep("^A.[.]R(c|n)", statenames), grep("^Ti[.]R(c|n)", statenames)), "novelprev"] <- 1
            
      tally[c(grep("^A.[.]R(rc|rn)", statenames), grep("^Ti[.]R(rc|cn)", statenames)), "rnovelprev"] <- 1
            
      tally[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames)),"inc"] <- #doens't include relapses (mostly early), but does include reinfections (mostly late) - i.3. this is a count of new infecitons, whether or not they are recognized as such
        apply(squaremat[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames)),  grep("^A",statenames)], 1, sum)
      
      tally[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames)),"rrinc"] <- #also only counts new infections with rr strain
        apply(squaremat[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames)),  grep("^A.[.]Rr",statenames)], 1, sum)
      
      tally[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames), grep("^R",statenames)),"rronsets"] <- #this alternative also includes resistance relapses (including resistance acquisitions with relapse)
        apply(squaremat[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames),grep("^R",statenames)),  grep("^A.[.]Rr",statenames)], 1, sum)
      
      tally[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames), grep("^R",statenames)),"panronsets"] <-
        apply(squaremat[c(grep("^S",statenames), grep("^C",statenames),grep("^L",statenames),grep("^R",statenames)),  grep("^A.[.]R(rc|rcn|rn)",statenames)], 1, sum)
      
      tally[grep("^R",statenames),"relapses"] <- 
        apply(squaremat[grep("^R",statenames),  grep("^A",statenames)], 1, sum)
      
      tally[c(grep("^A", statenames), grep("^Ti", statenames)), "tbdeaths"] <- rep(tbmort, each=(length(c(grep("^A", statenames), grep("^Ti", statenames)))/2))
      
      tally[c(grep("^A.[.]Rr", statenames), grep("^Ti[.]Rr", statenames)), "rrdeaths"] <- rep(tbmort, each=(length(c(grep("^A.[.]Rr", statenames), grep("^Ti[.]Rr", statenames)))/2))
      
      tally[c(grep("^A.{1,10}Hp$", statenames), grep("^T.{1,10}Hp$", statenames)), "hivtoo"] <- 1 # will later divide by prev        
      
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
      
      tally[grep("^Ts", statenames), "rxtime_s"] <- Tperiod.lengths[1:length(grep("^Ts", statenames))] #doesn't include ineffective (Ti) months
      tally[grep("^Tr", statenames), "rxtime_r"] <- Tperiod.lengths[1:length(grep("^Tr", statenames))] #doesn't include ineffective (Ti) months
      tally[grep("^Tn", statenames), "rxtime_n"] <- Tperiod.lengths[1:length(grep("^Tn", statenames))] #doesn't include ineffective (Ti) months
    }
    
    tallied <- t(tally) %*% state ; names(tallied) <- outcomes; 
    coprev <- tallied["hivtoo"]/tallied["prev"]; rrfrac <- tallied["rrinc"]/tallied["inc"]; tallied <- c(tallied, coprev, rrfrac); names(tallied) <- c(outcomes, "coprev", "rrfrac") 
    
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



#advance forward a bit:
## stateplus includes outcomes tracking variables
advance <- function(state, t0, addedt=1, rvary, nvary, reportsteps=1, fullpars) 
{
  o <- ode(state, seq(t0,t0+addedt,by=0.1), dxdt, parms=fullpars, rvary=rvary, nvary=nvary, do.tally=TRUE, method="adams")
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

