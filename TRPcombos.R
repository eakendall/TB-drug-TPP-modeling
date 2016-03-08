# # There's interest in understanding how much different improvements of a ds regimen 
# # can have an impact alone, if the efficacy doesn't improve (assuming noninferiority efficacy target.)
# 
# # For some characteristics (resistance related), the minimal target is worse than the SOC regimen. 
# # And for others (exclusions, scalability, and companion resistance if DST is being done), the
# # variations describe the fraction of patients who get the novel rather than SOC regimen (so only 
# # matters if other characteristics are better than SOC.)
# 
# # For DS,  I'll start with the following baseline:
#  efficacy minimal (=SOC)
#  barrier intermediate (~rif emergence with SOC)
#  companion intermediate (5%, likely if any existing drugs are used)
#  duration minimal (~SOC)
#  tolerability minimal (~SOC)
#  exclusions intermediate, HIV="nonHIV" (5% excluded, regardless of HIV status)
#  scalability minimal (only 50% of eligible patients get it)
#  
#  # And I'll model improvements to optimal (or for efficacy, intermediate; and in the case of barrier, companion, and exclusions, also worsening to minimal),
#  # and I'll model scenarios where improvements in tolerability and duration occur alone and
#  # where improvement from minimal to optimal in either tolerability or duration increase scalability by 25%,
#  # (and improvement to maximal in both tolerability and duration increases scalability by 50% to optimal or 100%)
# 
#  so the scenarios modeled will be:
#   baseline (minimals = c("efficacy", "duration", "tolerability", "scalability"))
#   efficacyint (minimals = c("duration", "tolerability", "scalability"))
#   barrierlow (minimals = c("efficacy", "duration", "tolerability", "scalability", "barrier"))
#   barrierhigh (minimals = c("efficacy", "duration", "tolerability", "scalability"), optimals= "barrier")
#   companionlow (minimals = c("efficacy", "duration", "tolerability", "scalability", "companion"))
#   companionhigh (minimals = c("efficacy", "duration", "tolerability", "scalability"), optimals= "companion")
#   tolerabilityhighalone (minimals = c("efficacy", "duration", "scalability"), optimals="tolerability")
#   tolerabilityhighwithscaleup (minimals = c("efficacy", "duration"), optimals="tolerability")
#   durationshortalone (minimals = c("efficacy", "tolerability", "scalability"), optimals="duration")
#   durationshortwithscaleup (minimals = c("efficacy", "tolerability"), optimals="duration")
#   durationandtolerabilityalone (minimals = c("efficacy", "scalability"), optimals=c("duration", "tolerability"))
#   durationandtolerabilitywithscaleup  (minimals = c"efficacy", optimals=c("duration", "tolerability", "scalability"))
#   exclusionslow (minimals = c("efficacy", "duration", "tolerability", "scalability", "exclusions"))
#   exclusionshigh (minimals = c("efficacy", "duration", "tolerability", "scalability"), optimals= "exclusions")
#   scalabilityhigh (minimals = c("efficacy", "duration", "tolerability"), optimals= "scalability")
#   
# ,  all with HIV="nonHIV"
#  

taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1] 
ntasks <- as.numeric(commandArgs(trailingOnly=TRUE))[2] 
tname <- commandArgs(trailingOnly=TRUE)[3]
targetpt <- commandArgs(trailingOnly=TRUE)[4]
DST <- commandArgs(trailingOnly=TRUE)[5]

location<-"../scratch/"
tag <- "20160304p"; Nsims_ds <- 250

rDSTall <- ifelse(targetpt=="DS", TRUE, FALSE)
ilimits <- ceiling(seq(0,Nsims_ds, length=ntasks+1))
currenttag <- paste0(tname,"_",tag)
if(rDSTall==TRUE) currenttag <- paste0("rDSTall.",currenttag)

source("TPPmat.R")

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
values <- set.values()
mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- set.novelvalues()$elementnames

alldrout <- numeric(0)
i <- 1; while(file.exists(paste0(location,"DRcalibration_",currenttag,".",i,".csv")))
{alldrout <- rbind(alldrout, read.csv(paste0(location,"DRcalibration_",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0

currenttag <- paste0(currenttag,".",taskid)

drout <- alldrout[alldrout$ids %in% (ilimits[taskid]+1):ilimits[taskid+1],]
tolerance <- 1.5
drout <- drout[drout[,"rrinc"]/drout[,"inc"] > 1/tolerance*drout[,"targetdr"] & drout[,"rrinc"]/drout[,"inc"] < tolerance*drout[,"targetdr"], ]  #within 3fold if rr incident fraction target

if (targetpt=="DS") scenarios <- c("baseline", "efficacyint", "barrierlow", "barrierhigh", "companionlow", "companionhigh", 
               "tolerabilityhighalone", "tolerabilityhighwithscaleup","durationshortalone", "durationshortwithscaleup", 
               "durationandtolerabilityalone", "durationandtolerabilitywithscaleup", 
               "exclusionslow", "exclusionshigh", "scalabilityhigh")
if (targetpt=="DR") stop("DR scenarios not yet defined!")


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

    minimals==""; optimals==""
    if (s=="baseline") { minimals = c("efficacy", "duration", "tolerability", "scalability") }
    if (s=="efficacyint") { minimals = c("duration", "tolerability", "scalability") }
    if (s=="barrierlow") { minimals = c("efficacy", "duration", "tolerability", "scalability", "barrier") }
    if (s=="barrierhigh") { minimals = c("efficacy", "duration", "tolerability", "scalability"); optimals= "barrier"}
    if (s=="companionlow") { minimals = c("efficacy", "duration", "tolerability", "scalability", "companion") }
    if (s=="companionhigh") { minimals = c("efficacy", "duration", "tolerability", "scalability"); optimals= "companion" }
    if (s=="tolerabilityhighalone")  { minimals = c("efficacy", "duration", "scalability"); optimals="tolerability" }
    if (s=="tolerabilityhighwithscaleup") { minimals = c("efficacy", "duration"); optimals="tolerability" }
    if (s=="durationshortalone") { minimals = c("efficacy", "tolerability", "scalability"); optimals="duration" }
    if (s=="durationshortwithscaleup") { minimals = c("efficacy", "tolerability"); optimals="duration" }
    if (s=="durationandtolerabilityalone") {minimals = c("efficacy", "scalability"); optimals=c("duration", "tolerability") }
    if (s=="durationandtolerabilitywithscaleup")  { minimals = "efficacy"; optimals=c("duration", "tolerability", "scalability") }
    if (s=="exclusionslow") { minimals = c("efficacy", "duration", "tolerability", "scalability", "exclusions") }
    if (s=="exclusionshigh") { minimals = c("efficacy", "duration", "tolerability", "scalability"); optimals= "exclusions" }
    if (s=="scalabilityhigh") { minimals = c("efficacy", "duration", "tolerability"); optimals= "scalability" }
    
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
