# want to look (for DS DSTall rDSTall for now) at:
# non-HIV exclusions (11% vs 5% vs none) with minimal, intermediate, and optimal efficacy, and with HIV exclusions = 0 
# HIV exclusions (all, 5%, none) with minimal, intermediate, and optimal efficacy, with non-HIV exclusions = 0
# 
# so 2 3x3 arrays, one for HIV and one for non-HIV exclusions

taskid <- 1#as.numeric(commandArgs(trailingOnly=TRUE))[1] #the idr's we want to run
tname <- "SouthAfrica"#commandArgs(trailingOnly=TRUE)[2]
targetpt <- "DR"#commandArgs(trailingOnly=TRUE)[3]
DST <- "DSTall"#commandArgs(trailingOnly=TRUE)[4]
rDSTall <- FALSE #commandArgs(trailingOnly=TRUE)[5]
location<-""#"../source/"

tag <- "20160111"
currenttag <- paste0(tname,"_",tag)
if (rDSTall == TRUE) currenttag <- paste0("rDSTall.",currenttag)
tasktag <- paste0(currenttag,".idr",taskid)
source("TPPmat.R")

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
values <- set.values()
genericvalues <- mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- c("all",set.novelvalues()$elementnames)

alldrout <- numeric(0)
i <- 1; while(file.exists(paste0(location,"DRcalibration_",currenttag,".",i,".csv")))
{alldrout <- rbind(alldrout, read.csv(paste0(location,"DRcalibration_",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0

drout <- alldrout[alldrout$idr==taskid,]
tolerance <- 1.5
drout <- drout[drout[,"rrinc"]/drout[,"inc"] > 1/tolerance*drout[,"targetdr"] & drout[,"rrinc"]/drout[,"inc"] < tolerance*drout[,"targetdr"], ] 


header <- c("inew", "ids","idr","targetprev","targetcoprev", "targetdr", "targetpt","DST", "rDSTall", names(unlist(genericvalues)))
header <- append(header, paste0( rep(tallynames, times=2*3*3), "10",
                         rep(c("HIV","nonHIV") ,each=3*3*length(tallynames)),
                         "exclusions", rep(c("minimal","intermediate","optimal"), each=3*length(tallynames)),
                         "efficacy", rep(c("minimal","intermediate","optimal"), each=length(tallynames)) ) )
                                                                 
if(!file.exists(paste0(location,"Exclusions","_", targetpt,DST,"_",tasktag,".csv"))) { write(header,  file=paste0(location,"Exclusions","_", targetpt,DST,"_",tasktag,".csv"), sep=",", ncol=length(header)) }
  
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

  for (H in c("HIV","nonHIV")) for (exclusions in c("minimal","intermediate","optimal")) for (efficacy in c("minimal","intermediate","optimal"))
  {
    print(paste0("Evaluating TRP variation HIV ",H,", exclusions ",exclusions,", efficacy ",efficacy," for Simulation #", inew," of ",nrow(drout)," for ",targetpt,DST,currenttag))
    valueset <- sampleTRP(mergedvalues = genericvalues, targetpt = targetpt, DST = DST, 
                          minimals=c("exclusions","efficacy")[which(c(exclusions,efficacy)=="minimal")], 
                          optimals=c("exclusions","efficacy")[which(c(exclusions,efficacy)=="optimal")],
                          HIV=H) # HIV = both (default), HIV, or nonHIV
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
  write(unlist(c(iter, valuevect, iresult)), file=paste0("Exclusions","_", targetpt,DST,"_",tasktag,".csv"), sep=",", append=TRUE, ncol=length(header))
}
