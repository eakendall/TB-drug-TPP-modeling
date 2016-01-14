source("TPPmat.R")
location <- ""
dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
values <- set.values()
genericvalues <- mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
tallynames <- readRDS(paste0(location,"tallynames_20160111.RDS"))
elementnames <- c("all", set.novelvalues()$elementnames)

wideheader <- c("inew", "ids","idr","targetprev","targetcoprev", "targetdr", "targetpt","DST", names(unlist(genericvalues)))
wideheader <- append(wideheader, paste0(rep(tallynames,times=11*3),rep(rep(0:10, each=length(tallynames)), times=3), 
                                        rep(c("allminimal", "allintermediate","alloptimal"), each=11*length(tallynames))))
for (i in 2:length(elementnames)) wideheader <- append(wideheader, 
                                                       paste0( rep(tallynames, times=2*11), 
                                                               rep( rep(0:10, each=length(tallynames)), 2),
                                                               rep(elementnames[i], each=22*length(tallynames) ),
                                                               rep( c("minimal","optimal"), each=11*length(tallynames) ) ) )

for (i in 1:25)
{
  r <- read.csv(paste0(location,"unfixed/TRPwideoutput_DSDSTnone_rDSTall.India_20160111.",i,".csv"))
  r2 <- cbind(r[seq(2,nrow(r),by=2),], r[1+seq(2,nrow(r),by=2),])
  ncol(r2)
  r2 <- r2[,1:(ncol(r)+2)]
  head(r2)[,(-100:0) + ncol(r2)]
  r2[,which(colnames(r2)=="prev0allminimal"):(ncol(r2)-2)] <- 
    r2[,which(colnames(r2)=="cprev0allminimal"):(ncol(r2))]
  r2 <- r2[,1:(ncol(r2)-2)]
  r <- rbind(r[1,1:ncol(r2)],r2)
  
  
  colnames(r) <- wideheader
  write(c(wideheader, t(r)), sep=",", ncol=length(wideheader), file=paste0(location,"fixed/TRPwideoutput_DSDSTnone_rDSTall.India_20160111.",i,".csv"))
}