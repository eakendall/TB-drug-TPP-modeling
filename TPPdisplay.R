source("TPPmat.R")

currenttag <- "India_20160105"; tolerance <- 1.5; location=""
# currenttag <- "Philippines_20151227"
# targetepi <- "Philippines"
# tolerance <- 2 #1.5 for India and SouthAfrica, 2 for Brazil and Philippines
# location="fromMARCC/Phil20/"

levels <- c("minimal","intermediate","optimal"); 
elementnames <- c("all", set.novelvalues()$elementnames)

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
grays <- c("gray30","gray60","gray90")
values <- set.values(); genericvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))

# wideheader <- c("inew", "ids","idr","targetprev","targetcoprev", "targetdr", "targetpt","DST")
# wideheader <- append(wideheader, paste0(rep(tallynames,times=11*3),rep(rep(0:10, each=length(tallynames)), times=3), 
#                                         rep(c("allminimal", "allintermediate","alloptimal"), each=11*length(tallynames))))
# for (i in 2:length(elementnames)) wideheader <- append(wideheader, 
#                                                        paste0( rep(tallynames, times=2*11), 
#                                                                rep( rep(0:10, each=length(tallynames)), 2),
#                                                                rep(elementnames[i], each=22*length(tallynames) ),
#                                                                rep( c("minimal","optimal"), each=11*length(tallynames) ) ) )
elementlabels <- c("All elements\nvaried", "% durably cured\n(optimal conditions)", "Regimen duration", 
                   "Prevalence of\nexisting resistance\nto regimen", "Barrier to\nacquired novel\ndrug resistance", 
                   "Exclusions and\ncontraindications", "Adherence/burden\nto patient")
dslabels <- c("","","", "94% cured","97% cured", "99% cured", "6 months","4 months", "2 months", 
              "10% resistant", "3% resistant", "None resistant", "1 acquired resistance per 20 treatments", "1 acquired resistance per 125 treatments", "Minimal acquired resistance",
              "Excludes 100% HIV; 11% non-HIV", "Excludes 10% HIV; 5% non-HIV", "No exclusions", "Same as standard of care", "1.5% fewer dropouts", "3% fewer dropouts")
drlabels <- c("","","", "76% cured","94% cured", "97% cured", "20 months","9 months", "6 months", 
              "15% resistant", "5% resistant", "None resistant", "1 acquired resistance per 10 treatments", "1 acquired resistance per 20 treatments", "1 acquired resistance per 125 treatments",
              "Excludes 100% HIV; 11% non-HIV", "Excludes 10% HIV; 5% non-HIV", "No exclusions", "Same as standard of care", "3% fewer dropouts", "6% fewer dropouts")

####### pull in data ###########

drout <- numeric(0)
for (i in 1:10) drout <- rbind(drout, read.csv(paste0(location,"DRcalibration_",currenttag,".",i,".csv"), header = TRUE)) #saved results from dr sampling runs at time 0
drout <- drout[ drout[,"rrinc"]/drout[,"inc"] > 1/tolerance*drout[,"targetdr"] & drout[,"rrinc"]/drout[,"inc"] < tolerance*drout[,"targetdr"], ]  #within 3fold if rr incident fraction target


allnovelwide <- list()
for (targetpt in c("DS","DR")) for (DST in c("DSTall","DSTnone")) #for (targetepi in names(targetepis))
{
  i <- 1; allnovelwide[[paste0(targetpt, DST)]] <- numeric(0)
  while (file.exists(paste0(location,"TRPwideoutput_", targetpt, DST,"_", currenttag,".",i,".csv")))
  {
    allnovelwide[[paste0(targetpt, DST)]] <- rbind(allnovelwide[[paste0(targetpt, DST)]], read.csv(paste0(location,"TRPwideoutput_", targetpt, DST,"_", currenttag,".",i,".csv")))
    i <- i+1
  }
}  


###############################
# Once country data loaded 

# need to specificy targetpt and DST here

novelwide <- allnovelwide[["DSDSTall"]]

outcome <- c("tbdeaths") #can set up loop over multiple outcomes


traj <- array(0, dim=c(11,3,5)); dimnames(traj) <- list("t"=0:10, "level"=levels, "q"=c(0.025,0.25,0.5,0.075,0.975))

for (t in 0:10) for (l in levels) traj[t+1,l,] <- quantile(novelwide[,colnames(novelwide)==paste0(outcome, t, "all",l)], c(0.025,0.25,0.5,0.075,0.975))

final_abs <- final_diff <- final_pct <- array(0,dim=c( length(elementnames) , 2 , 5 )); 
  final_down <- final_pctdown <- array(0,dim=c( length(elementnames) , 3 , 5 )); 
dimnames(final_abs) <- dimnames(final_diff) <- dimnames(final_pct) <-  list("vary"=elementnames, "level"=c("minimal", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))
  dimnames(final_down) <- dimnames(final_pctdown) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))

for (vary in elementnames) 
{ 
  final_abs[vary,1,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"minimal")], c(0.025,0.25,0.5,0.075,0.975))
  final_abs[vary,2,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"optimal")], c(0.025,0.25,0.5,0.075,0.975))
  
  final_diff[vary,1,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"minimal")] - novelwide[ , paste0(outcome,"10allintermediate")], c(0.025,0.25,0.5,0.075,0.975))
  final_diff[vary,2,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"optimal")] - novelwide[ , paste0(outcome,"10allintermediate")], c(0.025,0.25,0.5,0.075,0.975))
  
  final_pct[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - novelwide[ , paste0(outcome,"10allintermediate")] )/
                                    novelwide[ , paste0(outcome,"10allintermediate")], c(0.025,0.25,0.5,0.075,0.975))
  final_pct[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - novelwide[ , paste0(outcome,"10allintermediate")] )/
                                    novelwide[ , paste0(outcome,"10allintermediate")], c(0.025,0.25,0.5,0.075,0.975))

  final_down[vary,1,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_down[vary,2,] <- quantile(novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_down[vary,3,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  
  final_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
}  

par(mar=c(4,4,1,1), mfrow=c(1,1))
plot(0:10, traj[,"minimal","0.5"], ylim=c(0,targetepis[[targetepi]][1]/8), type='l', col="red", xlab="Years after introduction", ylab="Annual TB mortality")
points(0:10, traj[,"intermediate","0.5"], type='l', col='orange')
points(0:10, traj[,"optimal","0.5"], type='l', col='green')
points(c(0,10), c(median(drout[,"tbdeaths"]),median(drout[,"tbdeaths10"])), col='black', type='l') 
legend("bottomleft", legend= c("Continued current standard", "Minimal novel DS regimen", "Intermediate novel DS regimen", "Optimal novel DS regimen"),
       col=c("black","red","orange","green"), lty=1)


par(mar=c(3,5,4,3), mfrow=c(1,1))
bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="Percent reduction in year 10 TB mortality with novel regimen\n compared to projection under continued current practice", 
                    xlab="", las=2, cex.lab=1, 
                    ylim=c(120*final_pctdown[1,3,"0.5"], 8),
                    legend=c("minimal","intermediate","optimal"), args.legend=list(title="Level of varied element(s)",x="bottomright",cex=0.8),
                    main="", #paste0("Novel regimen for DS TB, with universal DST"), cex.main=1,
                    col=cols, names.arg=rep("", length(elementnames)))
text(bpctdown+0.4, -0.5, dslabels ,cex=0.9, pos=2, srt=90, col="black", font=2) 
text(colMeans(bpctdown)-0.5 ,1, elementlabels, cex=1, pos=4, srt=90, font=1)
mtext("Varied TRP element(s)", side=3, line=3)
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)

# Make this 2x2 for DS/DR and w/wo DST

novelwide <- allnovelwide[["DRDSTall"]]

outcome <- c("tbdeaths") #can set up loop over multiple outcomes
y <- 10

final_pctdown <- array(0,dim=c( length(elementnames) , 3 , 5 )); 
dimnames(final_pctdown) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))

for (vary in elementnames) 
{ 
  final_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, y, vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, y,"allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, y, vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
}  

par(mar=c(3,5,8,3), mfrow=c(1,2))
bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="Percent reduction in year 10 total TB mortality with novel regimen\n compared to projection under continued current practice", 
                    xlab="", las=2, cex.lab=1, ylim=c(-10,1),
                    legend=c("minimal","intermediate","optimal"), args.legend=list(title="Level of varied element(s)",x="bottomright",cex=0.8),
                    main="",
                    col=cols, names.arg=rep("", length(elementnames)))
text(bpctdown+0.4, -0.1, drlabels ,cex=0.8, pos=2, srt=90, col="black", font=2) 
text(colMeans(bpctdown) ,0.1, elementlabels, cex=1, pos=4, srt=90, xpd=NA)
mtext("Varied TRP element(s)", side=3, line=7)
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)

outcome <- c("rrdeaths") #can set up loop over multiple outcomes

final_pctdown <- array(0,dim=c( length(elementnames) , 3 , 5 )); 
dimnames(final_pctdown) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))

for (vary in elementnames) 
{ 
  final_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
}  

bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="Percent reduction in year 10 rifampin-resistant TB mortality", 
                    xlab="", las=2, cex.lab=1, ylim=c(-50,5),
                    legend=c("minimal","intermediate","optimal"), args.legend=list(title="Level of varied element(s)",x="bottomright",cex=0.8),
                    main="", 
                    col=cols, names.arg=rep("", length(elementnames)))
text(bpctdown+0.4, -0.1, drlabels ,cex=0.8, pos=2, srt=90, col="black", font=2) 
text(colMeans(bpctdown) ,0.8, elementlabels, cex=1, pos=4, srt=90, xpd=NA)
mtext("Varied TRP element(s)", side=3, line=7)
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)





############ duration

# sensitivity analysis: ltfu at 2 months
ltfu2mo <- rbind(read.csv("TRPwideoutput_DSDSTall_ltfu2mo.India_20160105.1.csv"), read.csv("TRPwideoutput_DSDSTall_ltfu2mo.India_20160105.2.csv"))
r <- read.csv("DRcalibration_ltfu2mo.India_20160105.1.csv", header=FALSE); s <- read.csv("DRcalibration_ltfu2mo.India_20160105.2.csv"); 
colnames(r) <- colnames(s); drltfu <- rbind(r,s)
plotpctdown(outcome="tbdeaths", scenario="DSDSTall", elements=elementnames, novelwide=ltfu2mo, drout =drltfu, barlabels=TRUE)
mtext("With all losses to follow up occurring at 2 months (to maximize the impact of regimen duration)", side=1, line=4, font=2)

# resource use: barplots as above but for outcomes of diagnoses, DSTs (rif and novel in same plot), and rxmonths (all 3 in same plot)

cumresource <- resource <- array(0,dim=c( length(elementnames) , 3 , 4 )); 
dimnames(cumresource) <- dimnames(resource) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "reg"=c("First-line","Second-line","Novel", "Treatment courses"))

for (vary in elementnames) for (nreg in 1:4) 
{ 
  outcome <- c("rxtime_s","rxtime_r","rxtime_n", "dxs")[nreg]
  resource[vary,1,nreg] <- median(novelwide[ , paste0(outcome, "10", vary,"minimal")])
  resource[vary,2,nreg] <- median(novelwide[ , paste0(outcome, "10allintermediate")])
  resource[vary,3,nreg] <- median(novelwide[ , paste0(outcome, "10", vary,"optimal")])
  for (t in 1:10)
  {
    cumresource[vary,1,nreg] <- cumresource[vary,1,nreg] + median(novelwide[ , paste0(outcome, t, vary,"minimal")])
    cumresource[vary,2,nreg] <- cumresource[vary,2,nreg] + median(novelwide[ , paste0(outcome, t,"allintermediate")])
    cumresource[vary,3,nreg] <- cumresource[vary,3,nreg] + median(novelwide[ , paste0(outcome, t, vary,"optimal")])
  }
}  

cols <- c("pink", "beige","palegreen")

par(mar=c(4,5,3,2), mfrow=c(2,2), oma=c(2,1,1,1)) 
bres <- barplot(12*cbind(aperm(resource,c(1,3,2))["efficacy",1:3,], aperm(resource,c(1,3,2))["duration",1:3,], aperm(resource,c(1,3,2))["companion",1:3,]), beside = FALSE, 
                      space=c(0.75,0.25,0.25), 
                    col=cols, ylab="Patient-months of treatment\nin year 10, by regimen")
legend(x = bres[1], y=max(12*rowSums(resource["efficacy",,1:3]))+10, xjust=0.2, yjust=0, fill=rev(cols),
       c("novel regimen","second-line regimen","first-line regimen"), xpd=NA)
text(bres[c(2,5,8)], -10*12, c("efficacy","duration","companion\nresistance"), cex=1, pos=1, xpd=NA)

bres <- barplot(12*cbind(aperm(cumresource,c(1,3,2))["efficacy",1:3,], aperm(cumresource,c(1,3,2))["duration",1:3,], aperm(cumresource,c(1,3,2))["companion",1:3,]), beside = FALSE, 
                space=c(0.75,0.25,0.25), 
                col=cols, ylab="Cumulative patient-months of treatment\nthrough year 10, by regimen")
text(bres[c(2,5,8)], -100*12, c("efficacy","duration","companion\nresistance"), cex=1, pos=1, xpd=NA)

bres <- barplot(c(aperm(resource,c(1,3,2))["efficacy",4,], aperm(resource,c(1,3,2))["duration",4,], aperm(resource,c(1,3,2))["companion",4,]), beside = TRUE, 
                space=c(0.75,0.25,0.25),
                col=c("gray30","gray60","gray90"), ylab="Total treatment courses initiated", xlab="")
text(bres[c(2,5,8)], -25, c("efficacy","duration","companion\nresistance"), cex=1, pos=1, xpd=NA)


bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,c(2,3,4),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality", 
                    xlab="", cex.lab=1, 
                    ylim=c(100*final_pctdown[1,3,"0.5"], 8))
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,c(2,3,4),"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,c(2,3,4),"0.975"], angle=90, code=3, length=0.05)



###############
#Companion drug resistance

novelwide <- allnovelwide[["DSDSTnone"]]; colnames(novelwide) <- wideheader
outcome <- c("tbdeaths") #can set up loop over multiple outcomes

final_pctdown2 <- array(0,dim=c( length(elementnames) , 3 , 5 ));
dimnames(final_pctdown2) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))

for (vary in elementnames) 
{ 
  final_pctdown2[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown2[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown2[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
}  


par(mar=c(4,4,2,1), mfrow=c(1,2), oma=c(1,1,1,1))

bpctdown <- barplot(height = 100*aperm(final_pctdown2, c(2,1,3))[,c(4),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality", 
                    xlab="Baseline prevalence of\ncompanion-drug resistance", names.arg=c("10%","3%","0%"), cex.lab=1, main="Without DST for\nnovel regimen",
                    ylim=c(-15,4), col=c("gray30","gray60","gray90"))
arrows(bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,4,"0.025"], bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,4,"0.975"], angle=90, code=3, length=0.05)


library("sensitivity")
pcc(X = drout[,38:70], y= (novelwide[,"tbdeaths10companionminimal"] - drout[,"tbdeaths10"] )/ drout[,"tbdeaths10"], rank=TRUE)

###???!!! what is going on with dxrate here?? 




novelwide <- allnovelwide[["DRDSTall"]]

outcome <- c("rrdeaths") #can set up loop over multiple outcomes

final_pctdown2 <- array(0,dim=c( length(elementnames) , 3 , 5 ));
dimnames(final_pctdown2) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))

for (vary in elementnames) 
{ 
  final_pctdown2[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown2[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown2[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
}  

novelwide <- allnovelwide[["DRDSTnone"]]

outcome <- c("rrdeaths") #can set up loop over multiple outcomes

final_pctdown3 <- array(0,dim=c( length(elementnames) , 3 , 5 ));
dimnames(final_pctdown3) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))

for (vary in elementnames) 
{ 
  final_pctdown3[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown3[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown3[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
}  


par(mfrow=c(1,2))
bpctdown <- barplot(height = 100*aperm(final_pctdown2, c(2,1,3))[,c(4),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality", 
                    xlab="Baseline prevalence of\ncompanion-drug resistance", names.arg=c("15%","5%","0%"), cex.lab=1,
                    col=c("gray30","gray60","gray90"), main="RR regimen\nWith novel regimen DST")
arrows(bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,4,"0.025"], bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,4,"0.975"], angle=90, code=3, length=0.05)


bpctdown <- barplot(height = 100*aperm(final_pctdown3, c(2,1,3))[,c(4),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality", 
                    xlab="Baseline prevalence of\ncompanion-drug resistance", names.arg=c("15%","5%","0%"), cex.lab=1,
                    col=c("gray30","gray60","gray90"), main="RR regimen\nWithout novel regimen DST")
arrows(bpctdown, aperm(100*final_pctdown3, c(2,1,3))[,4,"0.025"], bpctdown, aperm(100*final_pctdown3, c(2,1,3))[,4,"0.975"], angle=90, code=3, length=0.05)





###########################
## Barrier to resistance ##

outcome <- "panronsets"

novelwide <- allnovelwide[["DSDSTall"]]

final_pctdown <- array(0,dim=c( length(elementnames) , 3 , 5 ));
dimnames(final_pctdown) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))

for (vary in elementnames) 
{ 
  final_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
}  

novelwide <- allnovelwide[["DSDSTnone"]]; colnames(novelwide) <- wideheader

final_pctdown2 <- array(0,dim=c( length(elementnames) , 3 , 5 ));
dimnames(final_pctdown2) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))

for (vary in elementnames) 
{ 
  final_pctdown2[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown2[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown2[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
}  

par(mar=c(4,4,2,1), mfrow=c(1,2), oma=c(1,1,1,1))
bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    main="Novel regimen for DS TB,\nwith novel regimen DST",
                    col=c("gray30","gray60","gray90"))
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)

bpctdown <- barplot(height = 100*aperm(final_pctdown2, c(2,1,3))[,c(4),"0.5"], beside = TRUE, 
                    main="Novel regimen for DS TB,\nwithout novel regimen DST",
                    col=c("gray30","gray60","gray90"))
arrows(bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,4,"0.025"], bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,4,"0.975"], angle=90, code=3, length=0.05)




############
# exclusions

# synergy with maximal efficacy
synergyheader <- c("inew", "ids","idr","targetprev","targetcoprev", "targetdr", "targetpt","DST", names(unlist(genericvalues)))
for (i in (2:length(elementnames))) synergyheader <- append(synergyheader,  #fixed this after the fact
                                                            paste0( rep(tallynames, times=2*11), 
                                                                    rep( rep(0:10, each=length(tallynames)), 2),
                                                                    rep(elementnames[i], each=22*length(tallynames) ),
                                                                    rep( c("minimal","optimal"), each=11*length(tallynames) ) ) )

novelwide <- read.csv("efficacyhigh_DSDSTall_India_20160105.1.csv"); colnames(novelwide) <- synergyheader
drout <- alldrout[1:nrow(novelwide),]
final_pctdown <- array(0,dim=c( length(elementnames)-1 , 3, 5 ));
dimnames(final_pctdown) <- list("vary"=elementnames[-1], "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))

for (vary in elementnames[-1]) 
{ 
  final_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10efficacyoptimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  final_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
}  

par(mfrow=c(1,1))
bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    main=paste0("Novel ",substr(scenario,1,2)," TB regimen,\n",scenarionames[[scenario]]," DST"),
                    ylim=100*min(final_pctdown[,,"0.5"])* c(1.5,-0.2), ylab=outcomenames[[outcome]],
                    col=cols, xlab="")
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)
text(bpctdown+0.4, -0.1, dslabels[which(elements==elementnames) * 3 - (2:0)] ,cex=1, pos=2, srt=90, col="black", font=2)
mtext(paste0("Varying ",elementlabels[which(elements==elementnames)], "with fixed maximal efficacy"), side=1, cex=0.8, line=5)
}  
