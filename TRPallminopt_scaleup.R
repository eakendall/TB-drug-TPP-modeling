taskid <- as.numeric(commandArgs(trailingOnly=TRUE))[1]
ntasks <- as.numeric(commandArgs(trailingOnly=TRUE))[2]
tname <- commandArgs(trailingOnly=TRUE)[3]
targetpt <- commandArgs(trailingOnly=TRUE)[4]
DST <- commandArgs(trailingOnly=TRUE)[5]

location<-"../scratch/"
tag <- "20160419p"
currenttag <- paste0(tname,"_",tag)
if (targetpt=="DS") rDSTall <- TRUE else rDSTall <- FALSE #commandArgs(trailingOnly=TRUE)[5]
drtag <- ifelse(rDSTall == TRUE, paste0("rDSTall.",currenttag), currenttag)
tasktag <- paste0(currenttag,".",taskid)
source("TPPmat.R")

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
values <- set.values()
genericvalues <- mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- c("all",set.novelvalues()$elementnames)
if (targetpt=="DS") elementnames <- elementnames[1:8] else elementnames <- elementnames[c(1:7,9)] 

alldrout <- numeric(0)
i <- 1; while(file.exists(paste0(location,"DRcalibration_",drtag,".",i,".csv")))
{alldrout <- rbind(alldrout, read.csv(paste0(location,"DRcalibration_",drtag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0

tolerance <- 1.5
drout <- alldrout[alldrout[,"rrinc"]/alldrout[,"inc"] > 1/tolerance*alldrout[,"targetdr"] & alldrout[,"rrinc"]/alldrout[,"inc"] < tolerance*alldrout[,"targetdr"], ] 

ilimits <- ceiling(seq(0,nrow(drout), length=ntasks+1))

header <- c("inew", "ids","idr","targetprev","targetcoprev", "targetdr", "targetpt","DST", "rDSTall", names(unlist(genericvalues)))
header <- append(header, paste0( rep(tallynames, times=2*11), rep(0:10,each=length(tallynames)), rep(c("allmin", "allopt"), each=11*length(tallynames)) ) )

if(!file.exists(paste0(location,"Scaleupminopt","_", targetpt,DST,"_",tasktag,".csv"))) { write(header,  file=paste0(location,"Scaleupminopt","_", targetpt,DST,"_",tasktag,".csv"), sep=",", ncol=length(header)) }

for (inew in (ilimits[taskid]+1):ilimits[taskid+1])
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
  
  for (level in c("minimal","optimal"))
  {
    print(paste0("Evaluating TRP optimal for all ",level," for Simulation #", inew," of ",nrow(drout)," for ",targetpt,DST,tasktag))
    m <- o <- NA
    if (level=="minimal") {m <- elementnames[8]; o <- elementnames[2:7]}
    if (level=="optimal") {m <- NA; o <- elementnames[2:8]}
    valueset <- sampleTRP(mergedvalues = genericvalues, targetpt = targetpt, DST = DST, 
                                                   minimals=m, optimals=o , HIV="nonHIV")
    
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
    
    outset <- ode(y=unlist(novelstate), times=0:10, func=dxdt, parms=parset$fullpars, do.tally=TRUE, method="adams")
    
    iresult <- append(iresult, as.vector(t(outset[,tallynames])))
  }
  write(unlist(c(iter, valuevect, iresult)), file=paste0(location,"Scaleupminopt","_", targetpt,DST,"_",tasktag,".csv"), sep=",", append=TRUE, ncol=length(header))
}
