# Next steps:
#   
#   make sure WR is les for 2 vs 4 mo regimen ! 


sample <- read.csv("TRPwideoutput_DSDSTall_rDSTall.India_20160309p.1.csv")
samplestate <- read.csv("TRPstateoutput_DSDSTall_rDSTall.India_20160309p.1.csv")
sampledr <- read.csv("DRcalibration_rDSTall.India_20160309p.1.csv")
colnames(samplestate)

cbind(samplestate[samplestate$vary=="duration"&samplestate$level=="m"&samplestate$time==10,c("level","R.R0.Hn")], 
      samplestate[samplestate$vary=="all"&samplestate$level=="i"&samplestate$time==10,c("level","R.R0.Hn")],
      samplestate[samplestate$vary=="duration"&samplestate$level=="o"&samplestate$time==10,c("level","R.R0.Hn")])

# -- no, it's not!  Relapsers go back up ~5-10% as the regimen moves from 4 mo to 2 mo. Why?

cbind(samplestate[samplestate$vary=="all"&samplestate$level=="m"&samplestate$time==10,c("level","R.R0.Hn")], 
      samplestate[samplestate$vary=="all"&samplestate$level=="i"&samplestate$time==10,c("level","R.R0.Hn")],
      samplestate[samplestate$vary=="all"&samplestate$level=="o"&samplestate$time==10,c("level","R.R0.Hn")])
# differences here are huge and in the correct direction.

cbind(samplestate[samplestate$vary=="duration"&samplestate$level=="m"&samplestate$time==10,c("level","R.Rr.Hn")], 
      samplestate[samplestate$vary=="all"&samplestate$level=="i"&samplestate$time==10,c("level","R.Rr.Hn")],
      samplestate[samplestate$vary=="duration"&samplestate$level=="o"&samplestate$time==10,c("level","R.Rr.Hn")])
# there are also slightly (~0.1%) fewer rif-resistant relapsers for the 4mo ds regimen than the 2 (or 6) mo novel DS regimens

cbind(samplestate[samplestate$vary=="duration"&samplestate$level=="m"&samplestate$time==10,c("level","R.Rn.Hn")], 
      samplestate[samplestate$vary=="all"&samplestate$level=="i"&samplestate$time==10,c("level","R.Rn.Hn")],
      samplestate[samplestate$vary=="duration"&samplestate$level=="o"&samplestate$time==10,c("level","R.Rn.Hn")])
# and same for novel drug resistance


## running the first couple with a 3 mo regimen. 

taskid <- 1; ntasks <- 1; tname <- "India"; targetpt <- "DS"; DST <- "DSTall"; location <- ""
tag <- "20160309p"; Nsims_ds <- 250
rDSTall <- ifelse(targetpt=="DS", TRUE, FALSE)
ilimits <- ceiling(seq(0,Nsims_ds, length=ntasks+1))
currenttag <- "3mo"
#if(rDSTall==TRUE) currenttag <- paste0("rDSTall.",currenttag)
source("TPPmat.R")

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
values <- set.values()
mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- set.novelvalues()$elementnames

alldrout <- read.csv("DRcalibration_rDSTall.India_20160309p.1.csv")
drout <- alldrout[alldrout$ids %in% (ilimits[taskid]+1):ilimits[taskid+1],]
tolerance <- 1.5
drout <- drout[drout[,"rrinc"]/drout[,"inc"] > 1/tolerance*drout[,"targetdr"] & drout[,"rrinc"]/drout[,"inc"] < tolerance*drout[,"targetdr"], ]  #within 3fold if rr incident fraction target

set.novelvalues <- function()
{
  selections <- list()
  selections$poor_n <- array(c(0.06, 0.03, 0, 
                               0.24, 0.12, 0.06), dim=c(3,2))
  selections$months_n <- array(c(4,3,1, 
                                 20, 9, 6), dim=c(3,2))
  selections$cres <- array(c(0.1,0.1,  0.03,0.03,  0,0, 
                             0.15,0.15,  0.05,0.05,  0,0), dim=c(2,3,2))
  selections$barrierbase <- array(c(0.05, 0.008, 0,
                                    0.1, 0.05, 0.008), dim=c(3,2))
  selections$eligibility <- array(c(0.89,0,  0.95,0.9,  1,1,
                                    0.89,0,  0.95,0.9,  1,1), dim=c(2,3,2))
  
  
  selections$ltfurate_factor <- array(c(1.0,0.75,0.5,
                                        1.0,0.75,0.5), dim=c(3,2))
  
  selections$availability <- array(c(0.5,0.75,1,
                                     0.5,0.75,1), dim=c(3,2))
  selections$rifdx_increase <- c(0, 0.5, 1) # for DR only
  
  elementnames <- c("efficacy", "duration", "companion", "barrier", "exclusions", "tolerability", "uptake", "riftest")
  return(list("selections"=selections, "elementnames"=elementnames))
}

evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, targetpt=targetpt, DST=DST, tag=currenttag, rDSTall=rDSTall, location=location, ids=1) # can also specify ids and idr to run just a subset of drout

teststate <- read.csv("TRPstateoutput_DSDSTall_3mo.csv")
testwide  <- read.csv("TRPwideoutput_DSDSTall_3mo.csv")

samplestate <- samplestate[1:nrow(teststate), ]

cbind(samplestate[samplestate$vary=="duration"&samplestate$level=="m"&samplestate$time==10,c("level","R.R0.Hn")], #6mo
      samplestate[samplestate$vary=="all"&samplestate$level=="i"&samplestate$time==10,c("level","R.R0.Hn")], #4mo
      teststate[teststate$vary=="all"&teststate$level=="i"&teststate$time==10,c("level","R.R0.Hn")], #3mo
      samplestate[samplestate$vary=="duration"&samplestate$level=="o"&samplestate$time==10,c("level","R.R0.Hn")], #2mo
      teststate[teststate$vary=="duration"&teststate$level=="o"&teststate$time==10,c("level","R.R0.Hn")] #1mo
      )

# so the minimum relapsers occurs at 3 mo, and 1mo >> 2mo > 3mo
# m 6.217735     i 5.581920     i 5.117999     o 6.029117     o  7.481936


cbind(samplestate[samplestate$vary=="duration"&samplestate$level=="m"&samplestate$time==10,c("level","R.Rr.Hn")], #6mo
      samplestate[samplestate$vary=="all"&samplestate$level=="i"&samplestate$time==10,c("level","R.Rr.Hn")], #4mo
      teststate[teststate$vary=="all"&teststate$level=="i"&teststate$time==10,c("level","R.Rr.Hn")], #3mo
      samplestate[samplestate$vary=="duration"&samplestate$level=="o"&samplestate$time==10,c("level","R.Rr.Hn")], #2mo
      teststate[teststate$vary=="duration"&teststate$level=="o"&teststate$time==10,c("level","R.Rr.Hn")] #1mo
)

# differs  like 3.1606088     i 3.1565286     i 2.9802911     o 3.1604723     o 2.9925964

cbind(sample$tbdeaths10durationminimal, 
      sample$tbdeaths10allintermediate,
      testwide$tbdeaths10allintermediate, 
      sample$tbdeaths10durationoptimal,
      testwide$tbdeaths10durationoptimal)
# mortality ends up with roughly  the same pattern
# tbdeaths 28.54442 28.26316 27.83810 28.47360 28.91695

cbind(sample$inc10durationminimal, 
      sample$inc10allintermediate,
      testwide$inc10allintermediate, 
      sample$inc10durationoptimal,
      testwide$inc10durationoptimal)
# and inc not surprisingly differs by a small multiple of diff(relapse)
# 152.4468 151.5443 150.6285 152.2403 154.1707

cbind(samplestate[samplestate$vary=="duration"&samplestate$level=="m"&samplestate$time==10,c("level","Ts1.R0.Hn")], #6mo
      samplestate[samplestate$vary=="all"&samplestate$level=="i"&samplestate$time==10,c("level","Ts1.R0.Hn")], #4mo
      teststate[teststate$vary=="all"&teststate$level=="i"&teststate$time==10,c("level","Ts1.R0.Hn")], #3mo
      samplestate[samplestate$vary=="duration"&samplestate$level=="o"&samplestate$time==10,c("level","Ts1.R0.Hn")], #2mo
      teststate[teststate$vary=="duration"&teststate$level=="o"&teststate$time==10,c("level","Ts1.R0.Hn")] #1mo
)

# for Ts1 at time 10, 6>4, 3<<4, 2>4, 1>3, 1<<2
# 14652     m  3.116001     i  3.079002     i  2.037512     o  3.106385     o  2.132367

cbind(samplestate[samplestate$vary=="duration"&samplestate$level=="m"&samplestate$time==10,c("level","Ti.R0.Hn")], #6mo
      samplestate[samplestate$vary=="all"&samplestate$level=="i"&samplestate$time==10,c("level","Ti.R0.Hn")], #4mo
      teststate[teststate$vary=="all"&teststate$level=="i"&teststate$time==10,c("level","Ti.R0.Hn")], #3mo
      samplestate[samplestate$vary=="duration"&samplestate$level=="o"&samplestate$time==10,c("level","Ti.R0.Hn")], #2mo
      teststate[teststate$vary=="duration"&teststate$level=="o"&teststate$time==10,c("level","Ti.R0.Hn")] #1mo
)

# Ti same but less exaggerated pattern


cbind(samplestate[samplestate$vary=="duration"&samplestate$level=="m"&samplestate$time==1,c("level","Ts1.R0.Hn")], #6mo
      samplestate[samplestate$vary=="all"&samplestate$level=="i"&samplestate$time==1,c("level","Ts1.R0.Hn")], #4mo
      teststate[teststate$vary=="all"&teststate$level=="i"&teststate$time==1,c("level","Ts1.R0.Hn")], #3mo
      samplestate[samplestate$vary=="duration"&samplestate$level=="o"&samplestate$time==1,c("level","Ts1.R0.Hn")], #2mo
      teststate[teststate$vary=="duration"&teststate$level=="o"&teststate$time==1,c("level","Ts1.R0.Hn")] #1mo
)

# difference begins and is substantial even at time 1



# I THINK Ti's ARE STAYING THERE FOREVER UNLESS THEY DIE OR SELF CURE? 
plot(0:10, samplestate$Ti.R0.Hn[which(samplestate$vary=="duration"&samplestate$level=="o")[1:11]])
plot(0:10, samplestate$Ts1.R0.Hn[which(samplestate$vary=="duration"&samplestate$level=="o")[1:11]])
plot(0:10, samplestate$Tn1.R0.Hn[which(samplestate$vary=="duration"&samplestate$level=="o")[1:11]])
# Looks like we'd reached an equilibrium, with the total in T1 roughly equal to the number in each treatment month 
# (so if only 1/50 are failing, that's about 50 months ~ 4 years worth of failures in the box, which is consistent with staying there until death of self cure)
# and sure enough, switching to using novel regimen with half the failures cute the number in Ti at equilibrium by nearly 1/2. 

# sensitivty analysis of 2 vs 4 mo diff for different parameters like relapserate
# change 2 to 1 mo
# look at different levels of efficacy and adherence
