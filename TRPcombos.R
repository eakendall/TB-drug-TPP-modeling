# # There's interest in understanding how much different improvements of a ds regimen 
# # can have an impact alone, if the efficacy doesn't improve (assuming noninferiority efficacy target.)
# 
# # For some characteristics (resistance related), the minimal target is worse than the SOC regimen. 
# # And for others (exclusions, uptake, and companion resistance if DST is being done), the
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
#  uptake minimal (only 50% of eligible patients get it)
#  
#  # And I'll model improvements to optimal (or for efficacy, intermediate; and in the case of barrier, companion, and exclusions, also worsening to minimal),
#  # and I'll model scenarios where improvements in tolerability and duration occur alone and
#  # where improvement from minimal to optimal in either tolerability or duration increase uptake by 25%,
#  # (and improvement to maximal in both tolerability and duration increases uptake by 50% to optimal or 100%)
# 
#  so the scenarios modeled will be:
#   baseline (minimals = c("efficacy", "duration", "tolerability", "uptake"))
#   efficacyint (minimals = c("duration", "tolerability", "uptake"))
#   barrierlow (minimals = c("efficacy", "duration", "tolerability", "uptake", "barrier"))
#   barrierhigh (minimals = c("efficacy", "duration", "tolerability", "uptake"), optimals= "barrier")
#   companionlow (minimals = c("efficacy", "duration", "tolerability", "uptake", "companion"))
#   companionhigh (minimals = c("efficacy", "duration", "tolerability", "uptake"), optimals= "companion")
#   tolerabilityhighalone (minimals = c("efficacy", "duration", "uptake"), optimals="tolerability")
#   tolerabilityhighwithscaleup (minimals = c("efficacy", "duration"), optimals="tolerability")
#   durationshortalone (minimals = c("efficacy", "tolerability", "uptake"), optimals="duration")
#   durationshortwithscaleup (minimals = c("efficacy", "tolerability"), optimals="duration")
#   durationandtolerabilityalone (minimals = c("efficacy", "uptake"), optimals=c("duration", "tolerability"))
#   durationandtolerabilitywithscaleup  (minimals = c"efficacy", optimals=c("duration", "tolerability", "uptake"))
#   exclusionslow (minimals = c("efficacy", "duration", "tolerability", "uptake", "exclusions"))
#   exclusionshigh (minimals = c("efficacy", "duration", "tolerability", "uptake"), optimals= "exclusions")
#   uptakehigh (minimals = c("efficacy", "duration", "tolerability"), optimals= "uptake")
# ,  all with HIV="nonHIV"
#  
# # For DR, we want the baseline to instead include intermediate efficacy (now edited to be < standard 1stline regimen), intermediate duration (9 mo),
# intermediate uptake (won't make this our primary variable for increasing reach), but *** minimal riftest ***, and otherwise same.
#  # And I'll model improvements of each characteristic to optimal, 
#  # and in the case of barrier, companion, exclusions, efficacy, and duration (without change in scaleup), also worsening to minimal. 
#  # and I'll model scenarios where improvements in tolerability and duration occur alone and
#  # where improvement from minimal to optimal in either tolerability or duration reduce RR underdiagnosis by 50% (riftest to intermediate),
#  # (and improvement to optimal in both tolerability and duration increases riftest to optimal (100%))
# 
# so the DR scenarios modeled will be:
#   baseline (minimals = c("tolerability", "riftest"))
#   efficacymin (minimals = c("efficacy", "tolerability", "riftest"))
#   efficacyopt (minimals = c("tolerability", "riftest"), optimals="efficacy)
#   barrierlow (minimals = c("tolerability", "riftest", "barrier"))
#   barrierhigh (minimals = c("tolerability", "riftest"), optimals= "barrier")
#   companionlow (minimals = c("tolerability", "riftest", "companion"))
#   companionhigh (minimals = c("tolerability", "riftest"), optimals= "companion")
#   tolerabilityhighalone (minimals = c("riftest"), optimals="tolerability")
#   tolerabilityhighwithscaleup (minimals = "", optimals="tolerability")
#   durationlong (minimals = c("riftest", "duration"), optimals="")
#   durationshortalone (minimals = c("riftest"), optimals="duration")
#   durationshortwithscaleup (minimals = "", optimals="duration")
#   durationandtolerabilityalone (minimals = "riftest", optimals=c("duration", "tolerability"))
#   durationandtolerabilitywithscaleup  (minimals = "", optimals=c("duration", "tolerability", "riftest"))
#   exclusionslow (minimals = c("tolerability", "riftest", "exclusions"))
#   exclusionshigh (minimals = c("tolerability", "riftest"), optimals= "exclusions")
#   uptakehigh (minimals = c("tolerability"), optimals= "riftest")
#   
# ,  all with HIV="nonHIV"

taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1] 
ntasks <- as.numeric(commandArgs(trailingOnly=TRUE))[2] 
tname <- commandArgs(trailingOnly=TRUE)[3]
targetpt <- commandArgs(trailingOnly=TRUE)[4]
DST <- commandArgs(trailingOnly=TRUE)[5]

location<-"../scratch/"
tag <- "20160309p"; Nsims_ds <- 250

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

if (targetpt=="DS") scenarios <- c("baseline", "efficacyint", "barrierlow", "barrierhigh", "companionlow", "companionhigh", 
               "tolerabilityhighalone", "tolerabilityhighwithscaleup","durationshortalone", "durationshortwithscaleup", 
               "durationandtolerabilityalone", "durationandtolerabilitywithscaleup", 
               "exclusionslow", "exclusionshigh", "uptakehigh")
if (targetpt=="DR") scenarios <- c("baseline", "efficacymin", "efficacyopt", "barrierlow", "barrierhigh", "companionlow", "companionhigh", 
                                   "tolerabilityhighalone", "tolerabilityhighwithscaleup","durationlong", "durationshortalone", "durationshortwithscaleup", 
                                   "durationandtolerabilityalone", "durationandtolerabilitywithscaleup", 
                                   "exclusionslow", "exclusionshigh", "uptakehigh")

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

    minimals<-""; optimals<-""

    if (targetpt=="DS")
    {
      if (s=="baseline") { minimals = c("efficacy", "duration", "tolerability", "uptake") }
      if (s=="efficacyint") { minimals = c("duration", "tolerability", "uptake") }
      if (s=="barrierlow") { minimals = c("efficacy", "duration", "tolerability", "uptake", "barrier") }
      if (s=="barrierhigh") { minimals = c("efficacy", "duration", "tolerability", "uptake"); optimals= "barrier"}
      if (s=="companionlow") { minimals = c("efficacy", "duration", "tolerability", "uptake", "companion") }
      if (s=="companionhigh") { minimals = c("efficacy", "duration", "tolerability", "uptake"); optimals= "companion" }
      if (s=="tolerabilityhighalone")  { minimals = c("efficacy", "duration", "uptake"); optimals="tolerability" }
      if (s=="tolerabilityhighwithscaleup") { minimals = c("efficacy", "duration"); optimals="tolerability" }
      if (s=="durationshortalone") { minimals = c("efficacy", "tolerability", "uptake"); optimals="duration" }
      if (s=="durationshortwithscaleup") { minimals = c("efficacy", "tolerability"); optimals="duration" }
      if (s=="durationandtolerabilityalone") {minimals = c("efficacy", "uptake"); optimals=c("duration", "tolerability") }
      if (s=="durationandtolerabilitywithscaleup")  { minimals = "efficacy"; optimals=c("duration", "tolerability", "uptake") }
      if (s=="exclusionslow") { minimals = c("efficacy", "duration", "tolerability", "uptake", "exclusions") }
      if (s=="exclusionshigh") { minimals = c("efficacy", "duration", "tolerability", "uptake"); optimals= "exclusions" }
      if (s=="uptakehigh") { minimals = c("efficacy", "duration", "tolerability"); optimals= "uptake" }
    }
    if (targetpt=="DR")
    {
      if (s=="baseline") { minimals = c("tolerability", "riftest") }
      if (s=="efficacymin") { minimals = c("efficacy", "tolerability", "riftest") }
      if (s=="efficacyopt") { minimals = c("tolerability", "riftest"); optimals="efficacy" }
      if (s=="barrierlow") { minimals = c("tolerability", "riftest", "barrier") }
      if (s=="barrierhigh") { minimals = c("tolerability", "riftest"); optimals= "barrier" }
      if (s=="companionlow") { minimals = c("tolerability", "riftest", "companion") }
      if (s=="companionhigh") { minimals = c("tolerability", "riftest"); optimals= "companion" }
      if (s=="tolerabilityhighalone") { minimals = "riftest"; optimals="tolerability" }
      if (s=="tolerabilityhighwithscaleup") { minimals = ""; optimals="tolerability" }
      if (s=="durationlong") { minimals = c("riftest", "duration"); optimals="" }
      if (s=="durationshortalone") { minimals = "riftest"; optimals="duration" }
      if (s=="durationshortwithscaleup") { minimals = ""; optimals="duration" }
      if (s=="durationandtolerabilityalone") { minimals = "riftest"; optimals=c("duration", "tolerability") }
      if (s=="durationandtolerabilitywithscaleup")  { minimals = ""; optimals=c("duration", "tolerability", "riftest") }
      if (s=="exclusionslow") { minimals = c("tolerability", "riftest", "exclusions") }
      if (s=="exclusionshigh") { minimals = c("tolerability", "riftest"); optimals= "exclusions" }
      if (s=="uptakehigh") { minimals = "tolerability"; optimals= "riftest" }  
    }

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
