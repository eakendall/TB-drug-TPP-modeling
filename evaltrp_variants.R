#source the function below, then run this:

taskid <- 1#as.numeric(commandArgs(trailingOnly=TRUE))[1] #should have same array as dr runs
tname <- "India"#commandArgs(trailingOnly=TRUE)[2]
targetpt <- "DS"#commandArgs(trailingOnly=TRUE)[3]
DST <- "DSTall"#commandArgs(trailingOnly=TRUE)[4]
#startrow <- as.numeric(commandArgs(trailingOnly=TRUE))[5]
location<-""

tag <- "20160105"
currenttag <- paste0(tname,"_",tag,".",taskid)

source("TPPmat.R")

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
values <- set.values()
mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- set.novelvalues()$elementnames


drout <- screendrout(drout_filename=paste0("DRcalibration_",currenttag,".csv"), tolerance=1.5) # limits output to those simulations with DR fraction of incidence in targetepi range

#startrow <- max(read.csv(paste0(location,"TRPwideoutput_",targetpt,DST,"_",currenttag,".csv"))$inew) +1

evaltrp_synergy(synergy="efficacy",highlow="high",genericvalues = mergedvalues, drsetup = drsetup, drout=drout, rows=5:nrow(drout), targetpt=targetpt, DST=DST, tag=currenttag) # can also specify ids and idr to run just a subset of drout





evaltrp_synergy <- function(synergy="efficacy", highlow="high", genericvalues, drsetup, drout, ids, idr, rows, targetpt="DS", DST="DSTall", tag=currenttag,location="") # uses merged but not unlisted values
{
  if(missing(genericvalues)) {genericvalues <- readRDS(paste0("genericvalues_",tag,".RDS"))} # source of parameters that will have constant values (as saved at start of sampling)
  if(length(genericvalues[[1]][[1]]>1)) { genericvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]])) } #merge into single list if needed
  if(missing(drsetup)) {drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)}
  if(missing(drout)) {drout <- read.csv(file = paste0("DRcalibration_", tag, ".csv"), header = TRUE)} #includes i,i,targetepi,beta, hivrate, sampledpars(ds/dr), finalstate
  
  if(missing(rows))  if (missing(ids) || missing(idr)) {rows <- 1:nrow(drout)} else {rows <- (1:nrow(drout))[(drout[,"ids"] %in% ids) & (drout[,"idr"] %in% idr)]}
  novelsetup <- setup.model(DRera = TRUE, treatSL = TRUE, treatnovel = TRUE)
  elementnames <- c("all", set.novelvalues()$elementnames)
  
  wideheader <- c("inew", "ids","idr","targetprev","targetcoprev", "targetdr", "targetpt","DST", names(unlist(genericvalues)))
  for (i in (2:length(elementnames))-1) wideheader <- append(wideheader,  #fixed this after the fact
                                                         paste0( rep(tallynames, times=2*11), 
                                                                 rep( rep(0:10, each=length(tallynames)), 2),
                                                                 rep(elementnames[i], each=22*length(tallynames) ),
                                                                 rep( c("minimal","optimal"), each=11*length(tallynames) ) ) )
  
  
  if(!file.exists(paste0(location,synergy,highlow,"_", targetpt,DST,"_",tag,".csv"))) { write(wideheader,  file=paste0(location,synergy,highlow,"_", targetpt,DST,"_",tag,".csv"), sep=",", ncol=length(wideheader)) }
  
  for (inew in rows)
  {
    iter <- unlist(c(inew,unlist(drout[inew,c("ids", "idr", "targetprev","targetcoprev","targetdr")]), targetpt, DST)) #will include these labels as part of returned output
    
    valuevect <- drout[inew, 5+(1:length(unlist(genericvalues)))] ; names(valuevect) <- names(unlist(genericvalues))
    v <- 0; for (pname in names(genericvalues))
    { genericvalues[[pname]] <- unlist(valuevect[v+(1:length(genericvalues[[pname]]))]); #names(genericvalues[[pname]]) <- names()
      v <- v + length(genericvalues[[pname]]) # move forward to start of next par vector in sampled values
    }    
    
    state <- drout[inew,drsetup$statenames]
    newstate <- numeric(length(novelsetup$statenames)); names(newstate) <- novelsetup$statenames
    newstate[drsetup$statenames] <- unlist(state)
    
    iresult <- numeric(0) #will become full row of results for this dr row (inew), for all elements and levels
    for (vary in elementnames[-1])
    {
      print(paste0("Evaluating TRP variation in ",vary,  ",  Simulation #", inew," of ",nrow(drout)," for ",targetpt,DST,currenttag))
      valueset <- list()
      { if (highlow=="high" )
        {
          valueset$m <- sampleTRP(mergedvalues = genericvalues, targetpt = targetpt, DST = DST, minimals=vary, optimals=synergy[!(synergy %in% vary)])
          valueset$o <- sampleTRP(mergedvalues = genericvalues, targetpt = targetpt, DST = DST, optimals=c(vary, synergy[!(vary %in% synergy)]))
        }
        if (highlow=="low" )
        {
          valueset$m <- sampleTRP(mergedvalues = genericvalues, targetpt = targetpt, DST = DST, minimals=c(vary, synergy[!(vary %in% synergy)]))
          valueset$o <- sampleTRP(mergedvalues = genericvalues, targetpt = targetpt, DST = DST, optimals=vary, minimals=synergy[!(synergy %in% vary)])          
        }
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
        
        iresult <- append(iresult, as.vector(t(outset[,tallynames])))
      }
    }
    
    write(unlist(c(iter, valuevect, iresult)), file=paste0(location,synergy,highlow,"_", targetpt,DST,"_",tag,".csv"), sep=",", append=TRUE, ncol=length(wideheader))
    
  }
  return(0)
}
