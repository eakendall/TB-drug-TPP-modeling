# # There's interest in understanding how much different improvements of a ds regimen 
# # can have an impact alone, if the efficacy doesn't improve (assuming noninferiority efficacy target.)
# 
# We also want to know how impact changes with scaleup (DS, although this isn't really different from exlusions) or riftest (DR, a critical consideration).
# 
# # For some characteristics (resistance related), the minimal target is worse than the SOC regimen, so we'll instead use intermediate in our baseline (and do separate anaysis for variation in resistance and DST practices). 
# # And for others (exclusions, uptake, and companion resistance if DST is being done), the  variations describe the fraction of patients who
# # get the novel rather than SOC regimen (so only matter if other characteristics are better than SOC.)
# 

# same approach for DS and DR:
# for each of low/med/high efficacy: 
#   minimal duration, tolerability, and exclusions, with intermediate scaleup and (for DR) riftest
#   duration optimal/optimal with increased (optimal) scaleup/riftest
#   tolerability high/high with increased scaleup/riftest
#   non-HIV exclusions optimal/optimal with increased scaleup/riftest
#   (all with intermediate barrier and companion)

# (And will do separate resistance analysis)

taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1] 
ntasks <- as.numeric(commandArgs(trailingOnly=TRUE))[2] 
tname <- commandArgs(trailingOnly=TRUE)[3]
targetpt <- commandArgs(trailingOnly=TRUE)[4]
DST <- commandArgs(trailingOnly=TRUE)[5]

location<-"../scratch/"
tag <- "20160313p"; Nsims_ds <- 50

rDSTall <- ifelse(targetpt=="DS", TRUE, FALSE)
ilimits <- ceiling(seq(0,Nsims_ds, length=ntasks+1))
currenttag <- paste0(tname,"_",tag)
if(rDSTall==TRUE) currenttag <- paste0("rDSTall.",currenttag)

source("TPPmat.R")

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
values <- set.values()
genericvalues <- mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- set.novelvalues()$elementnames

alldrout <- numeric(0)
i <- 1; while(file.exists(paste0(location,"DRcalibration_",currenttag,".",i,".csv")))
{alldrout <- rbind(alldrout, read.csv(paste0(location,"DRcalibration_",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0

currenttag <- paste0(currenttag,".",taskid)

drout <- alldrout[alldrout$ids %in% (ilimits[taskid]+1):ilimits[taskid+1],]
tolerance <- 1.5
drout <- drout[drout[,"rrinc"]/drout[,"inc"] > 1/tolerance*drout[,"targetdr"] & drout[,"rrinc"]/drout[,"inc"] < tolerance*drout[,"targetdr"], ]  #within 3fold if rr incident fraction target

# for each of low/med/high efficacy: 
#   minimal duration, tolerability, and exclusions, with intermediate scaleup (uptake for DS, and riftest (and uptake fixed at int) for DR)
#   duration optimal/optimal with increased (optimal) uptake/riftest
#   tolerability high/high with increased uptake/riftest
#   non-HIV exclusions optimal/optimal with increased uptake/riftest
#   (all with intermediate barrier and companion)
scenarios <- c("em_dm_tm_xm_sm_","em_dm_tm_xm_si_","em_dm_tm_xm_so_", 
                 "em_do_tm_xm_sm_", "em_do_tm_xm_si_", "em_do_tm_xm_so_", 
                 "em_dm_to_xm_sm_", "em_dm_to_xm_si_", "em_dm_to_xm_so_", 
                 "em_dm_tm_xo_sm_", "em_dm_tm_xo_si_", "em_dm_tm_xo_so_",
               "ei_dm_tm_xm_sm_", "ei_dm_tm_xm_si_", "ei_dm_tm_xm_so_", 
                 "ei_do_tm_xm_sm_","ei_do_tm_xm_si_", "ei_do_tm_xm_so_", 
                 "ei_dm_to_xm_sm_", "ei_dm_to_xm_si_", "ei_dm_to_xm_so_", 
                 "ei_dm_tm_xo_sm_","ei_dm_tm_xo_si_", "ei_dm_tm_xo_so_",
               "eo_dm_tm_xm_sm_","eo_dm_tm_xm_si_", "eo_dm_tm_xm_so_",
                 "eo_do_tm_xm_sm_","eo_do_tm_xm_si_", "eo_do_tm_xm_so_", 
                 "eo_dm_to_xm_sm_","eo_dm_to_xm_si_", "eo_dm_to_xm_so_", 
                  "eo_dm_tm_xo_sm_","eo_dm_tm_xo_si_", "eo_dm_tm_xo_so_")

header <- c("inew", "ids","idr","targetprev","targetcoprev", "targetdr", "targetpt","DST", "rDSTall", names(unlist(genericvalues)))
header <- append(header, paste0( rep(tallynames, times=length(scenarios)), "10",
                                 rep(scenarios, each=length(tallynames)) ) )

if(!file.exists(paste0(location,"TRPcombos","_", targetpt,DST,"_",currenttag,".csv"))) { write(header,  file=paste0(location,"TRPcombos","_", targetpt,DST,"_",currenttag,".csv"), sep=",", ncol=length(header)) }


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
  
  for (s in scenarios)
  {
    print(paste0("Evaluating TRP scenario ",s," for Simulation #", inew," of ",nrow(drout)," for ",targetpt,DST,currenttag))

    minimals<-character(); optimals<-character()

    if (grepl("em_",s)) minimals <- c(minimals, "efficacy")
    if (grepl("eo_",s)) optimals <- c(optimals, "efficacy")
    if (grepl("dm_",s)) minimals <- c(minimals, "duration")
    if (grepl("do_",s)) optimals <- c(optimals, "duration")
    if (grepl("tm_",s)) minimals <- c(minimals, "tolerability")
    if (grepl("to_",s)) optimals <- c(optimals, "tolerability")
    if (grepl("xm_",s)) minimals <- c(minimals, "exclusions")
    if (grepl("xo_",s)) optimals <- c(optimals, "exclusions")
    if (grepl("sm_",s) & targetpt=="DS") minimals <- c(minimals, "uptake")
    if (grepl("so_",s) & targetpt=="DS") optimals <- c(optimals, "uptake")
    if (grepl("sm_",s) & targetpt=="DR") minimals <- c(minimals, "riftest")
    if (grepl("so_",s) & targetpt=="DR") optimals <- c(optimals, "riftest")
    
    valueset <- sampleTRP(mergedvalues = genericvalues, targetpt = targetpt, DST = DST, 
                          minimals=minimals, optimals=optimals, HIV="nonHIV")
  
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
    
    outset <- ode(y=unlist(novelstate), times=0:10, func=dxdt, parms=parset$fullpars, do.tally=TRUE, method="adams")[11,]
    
    iresult <- append(iresult, outset[tallynames])
    
  }

  write(unlist(c(iter, valuevect, iresult)), file=paste0(location,"TRPcombos","_", targetpt,DST,"_",currenttag,".csv"), sep=",", append=TRUE, ncol=length(header))
}
