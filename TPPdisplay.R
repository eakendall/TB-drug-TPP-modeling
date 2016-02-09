source("TPPmat.R")

currenttag <- "India_20160201"; tolerance <- 1.5; location=""
targetepi <- "India"

source("displayfunction.R")


###############################
par(mar=c(3,6,7,1), mfrow=c(1,1),oma=c(0,0,0,0))
bpctdown <- barplot(height = array(-c(1,3,5,2,3,4), dim=c(3,2)), beside = TRUE, 
                    ylab="Impact of novel regimen\n(mortality reduction)", ylim=c(-6,0), yaxt='n',
                    xlab="", las=1, cex.lab=1.5, 
                    legend=c("minimal","intermediate","optimal"), args.legend=list(title="Level of varied element(s)", x=8, y= 2.2, cex=1.2),
                    col=cols, names.arg=c("Over full range\n of possible novel regimens","Varying only\nspecified characteristic"), cex.names=1.1)
arrows(bpctdown[1,], -5.5, bpctdown[3,], -5.5, angle=90, code=3, length=0.1)


# Results overview #

novelwide <- allnovelwide[["DSDSTall_rDSTall."]]; drout <- droutds <- alldrDST[1:nrow(novelwide),]

outcome <- c("tbdeaths") #can set up loop over multiple outcomes

traj <- array(0, dim=c(11,3,5)); dimnames(traj) <- list("t"=0:10, "level"=levels, "q"=c(0.025,0.25,0.5,0.75,0.975))
for (t in 0:10) for (l in levels) traj[t+1,l,] <- quantile(novelwide[,colnames(novelwide)==paste0(outcome, t, "all",l)], c(0.025,0.25,0.5,0.75,0.975))
# par(mar=c(4,4,1,1), mfrow=c(2,1))
plot(0:10, traj[,"minimal","0.5"], ylim=c(0,targetepis[[targetepi]][1]/4), type='l', col="deeppink", xlab="Years after introduction", ylab="Annual TB mortality",lwd=2)
points(0:10, traj[,"intermediate","0.5"], type='l', col='orange',lwd=2)
points(0:10, traj[,"optimal","0.5"], type='l', col='green',lwd=2)
points(c(0,10), c(median(drout[,"tbdeaths"]),median(drout[,"tbdeaths10"])), col='black', type='l',lwd=2) 
# points(c(0,10), c(median(drout[,"tbdeaths"]),median(drout[,"tbdeaths10"])+1.5), col='black', type='l',lwd=2) 
legend("bottomleft", legend= c("Continued current standard", "Minimal novel DS regimen", "Intermediate novel DS regimen", "Optimal novel DS regimen"),
       col=c("black","deeppink","orange","green"), lty=1,lwd=2)

final_pctdown <- array(0,dim=c( length(elementnames) , 3 , 5 )); 
dimnames(final_pctdown) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in elementnames) 
{ final_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  final_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  final_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}  

par(mar=c(2,4,8,1), mfrow=c(1,1), oma=c(0,0,0,0))
bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="% reduction (median [IQR])",
                    xlab="", las=2, cex.lab=1, 
                    ylim=c(145*final_pctdown[1,3,"0.5"], 2), yaxt='n',
                    legend=c("minimal","intermediate","optimal"),
                    args.legend=list(title="Level of varied element(s)", x="bottom",cex=0.9),
                    space=c(0,0,0,1.5,0,0,rep(c(0.5,0,0),5)),
                    col=cols, names.arg=rep("", length(elementnames)))
axis(2, at=seq(-25,0,by=5), labels=paste0(seq(-25,0,by=5),"%"),las=2,cex.axis=0.8)
mtext("Reduction in year 10 TB mortality with a novel DS-TB regimen,\n compared to projection under current standard of care", cex=1.2, line=5, side=3)
mtext("Varied TRP element(s)", side=3, line=3)
# text((bpctdown+0.4)[,2], -0.5, dslabels[4:6] ,cex=0.9, pos=2, srt=90, col="black", font=2) 
text(colMeans(bpctdown) ,0.5, elementlabels, cex=0.8, pos=3, srt=0, font=1, xpd=NA)
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.25"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.75"], angle=90, code=3, length=0.05, xpd=NA)
# arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)
# calculate contribution of each element to total novel regimen impact

# ffracs <- (final_pctdown[2:7,3,3]-final_pctdown[2:7,1,3] )/sum(final_pctdown[2:7,3,3]-final_pctdown[2:7,1,3] )
# ffracs <- (final_pctdown[2:7,3,3]-final_pctdown[2:7,1,3] )/(final_pctdown[1,3,3]-final_pctdown[1,1,3] )
# fp <- numeric(6); for (i in 1:6) fp[i] <- sum(ffracs[1:i])
# fp2 <- c(0,fp[-6]) + (fp - c(0,fp[-6]))/2
# 
# par(mar=c(1,1,3,1)) 
# b <- barplot(array(ffracs, dim=c(6,1)), horiz = TRUE, beside=FALSE, las=2, font=2, xaxt='n',
#              legend.text=shortelementlabels[-1], col=rainbow(6), args.legend=list(x=0.2,y=1.6))
# mtext("Fraction of total variation in novel regimen's impact attributable to\nvariation in specified characteristic", side=3, font=2, cex=1.2)
# text(fp2,b+0.1*(rep(c(0,1),3)),paste0(round(100*ffracs,0),"%"), font=2)
# 

allbut_ds <- read.csv("Allbut_DSDSTall_India_20160201.csv")


contribs <- array(0,dim=c(7,5)); dimnames(contribs) <- list("element"=elementnames, "q"=c(0.025,0.25,0.5,0.75,0.975))

for (vary in elementnames[2:7])
{
  contribs[vary,] <- quantile( 
      ( novelwide[ , paste0(outcome, "10alloptimal")] - drout[ , paste0(outcome,"10")]  - 
        (allbut_ds[ , paste0(outcome, "10allbut",vary)] - drout[ , paste0(outcome,"10")]) ) /
          (novelwide[ , paste0(outcome, "10alloptimal")] - drout[ , paste0(outcome,"10")]) , c(0.025,0.25,0.5,0.75,0.975) )
}
contribs["all",] <- quantile( 
  ( novelwide[ , paste0(outcome, "10alloptimal")] - drout[ , paste0(outcome,"10")]  - 
      (novelwide[ , paste0(outcome, "10allminimal")] - drout[ , paste0(outcome,"10")] ) )/
    (novelwide[ , paste0(outcome, "10alloptimal")] - drout[ , paste0(outcome,"10")]) , c(0.025,0.25,0.5,0.75,0.975) )

par(mar=c(1,1,3,1), mfrow=c(1,1)) 
b <- barplot(contribs[,3], horiz = TRUE, beside=TRUE, las=2, font=2, xaxt='n',
             legend.text=shortelementlabels, col=rainbow(7), args.legend=list(x=1.5, y=8, cex=0.9))
mtext("Median fraction of the total mortality impact of an all-optimal novel regimen that is lost\nwhen one characteristic [or all characteristics] is reduced to its minimal value", side=3, line=0, font=2, cex=1.2, xpd=NA)
text(contribs[,3]+0.05, b, paste0(round(contribs[,3]*100, 1),"%"))




# Same for DR regimen
## resize plot area to wide ##
novelwide3 <- novelwide <- allnovelwide[["DRDSTall_"]]; drout <- droutdr <- alldrout[1:nrow(novelwide),]

outcome <- c("tbdeaths") 
r_pctdown <- array(0,dim=c( length(elementnames) , 3 , 5 )); y <- 10
dimnames(r_pctdown) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in elementnames) 
{ r_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, y, vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  r_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, y,"allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  r_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, y, vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}  

outcome <- c("rrdeaths"); y <- 10
rr_pctdown <- array(0,dim=c( length(elementnames) , 3 , 5 )); 
dimnames(rr_pctdown) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in elementnames) 
{ rr_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, y, vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  rr_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, y,"allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  rr_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, y, vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}  


par(mar=c(2,4,2,1), mfrow=c(1,2), oma=c(0,0,2,0))
ylow=300*r_pctdown[1,3,"0.5"]
bpctdown <- barplot(height = 100*aperm(r_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="% reduction, total TB mortality (median [IQR])", 
                    main="", xlab="", las=2, cex.lab=1, 
                    ylim=c(1,0)*ylow, yaxt='n',
                    legend=c("minimal","intermediate","optimal"), args.legend=list(title="Level of varied element(s)", x=bpctdown[2,7], y=-3,cex=0.8),
                    space=c(0,0,0,1.5,0,0,rep(c(0.5,0,0),5)),
                    col=cols, names.arg=rep("", length(elementnames)))
axis(2, at=c(-5,0), labels=paste0(c(5,0),"%"),las=2,cex.axis=0.8)
# text((bpctdown+0.4)[,2], -0.1, drlabels[4:6] ,cex=0.9, pos=2, srt=90, col="black", font=2) 
text(colMeans(bpctdown) ,ylow, shortelementlabels, cex=0.8, pos=4, srt=90, font=2, xpd=NA)
mtext("Varied TRP element(s)", side=1, line=1, xpd=NA, cex=0.8)
arrows(bpctdown, aperm(100*r_pctdown, c(2,1,3))[,,"0.25"], bpctdown, aperm(100*r_pctdown, c(2,1,3))[,,"0.75"], angle=90, code=3, length=0.05, xpd=NA)

ylow <- 130*rr_pctdown[1,3,"0.5"]
bpctdown <- barplot(height = 100*aperm(rr_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="% reduction, DR-TB mortality (median [IQR])", 
                    main="", xlab="", las=2, cex.lab=1, 
                    ylim=c(1,0)*ylow, yaxt='n',
                    space=c(0,0,0,1.5,0,0,rep(c(0.5,0,0),5)),
                    col=cols, names.arg=rep("", length(elementnames)))
axis(2, at=-seq(0,-ylow,by=5), labels=paste0(seq(0,-ylow,by=5),"%"),las=2,cex.axis=0.8)
# text((bpctdown+0.5)[,2], ylow/100, drlabels[4:6] ,cex=0.9, pos=2, srt=90, col="black", font=2) 
text(colMeans(bpctdown) ,ylow, shortelementlabels, cex=0.8, pos=4, srt=90, font=2, xpd=NA)
mtext("Varied TRP element(s)", side=1, line=1, cex=0.8)
arrows(bpctdown, aperm(100*rr_pctdown, c(2,1,3))[,,"0.25"], bpctdown, aperm(100*rr_pctdown, c(2,1,3))[,,"0.75"], angle=90, code=3, length=0.05)

mtext("Reduction in year 10 TB mortality with novel DR-TB regimen,\n compared to projection under current standard of care",side=3,outer=TRUE, cex=1.4, line=-1, xpd=NA) 




allbut_dr <- read.csv("Allbut_DRDSTall_India_20160201.csv")

contribs <- array(0,dim=c(7,5)); dimnames(contribs) <- list("element"=elementnames, "q"=c(0.025,0.25,0.5,0.75,0.975))

for (vary in elementnames[2:7])
{
  contribs[vary,] <- quantile( 
    ( novelwide[ , paste0(outcome, "10alloptimal")] - drout[ , paste0(outcome,"10")]  - 
        (allbut_ds[ , paste0(outcome, "10allbut",vary)] - drout[ , paste0(outcome,"10")]) ) /
      (novelwide[ , paste0(outcome, "10alloptimal")] - drout[ , paste0(outcome,"10")]) , c(0.025,0.25,0.5,0.75,0.975) )
}
contribs["all",] <- quantile( 
  ( novelwide[ , paste0(outcome, "10alloptimal")] - drout[ , paste0(outcome,"10")]  - 
      (novelwide[ , paste0(outcome, "10allminimal")] - drout[ , paste0(outcome,"10")] ) )/
    (novelwide[ , paste0(outcome, "10alloptimal")] - drout[ , paste0(outcome,"10")]) , c(0.025,0.25,0.5,0.75,0.975) )

par(mar=c(1,1,3,1), mfrow=c(1,1)) 
b <- barplot(contribs[,3], horiz = TRUE, beside=TRUE, las=2, font=2, xaxt='n',
             legend.text=shortelementlabels, col=rainbow(7), args.legend=list(x=1.5, y=8, cex=0.9))
mtext("Median fraction of the rifampin-resistant mortality impact of an all-optimal novel DR regimen that is lost\nwhen one characteristic [or all characteristics] is reduced to its minimal value", side=3, line=0, font=2, cex=1.2, xpd=NA)
text(contribs[,3]+0.05, b, paste0(round(contribs[,3]*100, 1),"%"))

# ffracs <- (rr_pctdown[2:7,3,3]-rr_pctdown[2:7,1,3] )/sum(rr_pctdown[2:7,3,3]-rr_pctdown[2:7,1,3] )
# fp <- numeric(6); for (i in 1:6) fp[i] <- sum(ffracs[1:i])
# fp2 <- c(0,fp[-6]) + (fp - c(0,fp[-6]))/2
# par(mar=c(1,1,3,1), mfrow=c(1,1), oma=c(0,0,0,0)) 
# b <- barplot(array(ffracs, dim=c(6,1)), horiz = TRUE, beside=FALSE, las=2, font=2, xaxt='n',
#              legend.text=shortelementlabels[-1], col=rainbow(6), args.legend=list(x=0.2,y=1.6))
# mtext("Fraction of total variation in novel regimen's MDR-TB mortality impact attributable to\nvariation in specified characteristic", side=3, font=2, cex=1.2)
# text(fp2,b+0.1*(rep(c(0,1),3)),paste0(round(100*ffracs,0),"%"), font=2)
# 
############ duration
# DS regimen:
# resource use: barplots as above but for outcomes of diagnoses, DSTs (rif and novel in same plot), and rxmonths (all 3 in same plot)
novelwide <- novelwide1 <- allnovelwide[["DSDSTall_rDSTall."]]; drout <- alldrDST[1:nrow(novelwide),]

cumresource <- array(0,dim=c( length(elementnames) , 3 , 3)); resource <- array(0,dim=c( length(elementnames) , 3 , 3,5)); 
dimnames(cumresource) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "reg"=c("First-line","Second-line","Novel")); dimnames(resource) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "reg"=c("First-line","Second-line","Novel"), "q"=c(0.025,0.25,0.5,0.75,0.975))

for (vary in elementnames) for (nreg in 1:3) 
{ 
  outcome <- c("rxtime_s","rxtime_r","rxtime_n", "dxs")[nreg]
  resource[vary,1,nreg,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"minimal")],c(0.025,0.25,0.5,0.75,0.975))
  resource[vary,2,nreg,] <- quantile(novelwide[ , paste0(outcome, "10allintermediate")],c(0.025,0.25,0.5,0.75,0.975))
  resource[vary,3,nreg,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"optimal")],c(0.025,0.25,0.5,0.75,0.975))
  for (t in 1:10)
  {
    cumresource[vary,1,nreg] <- cumresource[vary,1,nreg] + median(novelwide[ , paste0(outcome, t, vary,"minimal")])
    cumresource[vary,2,nreg] <- cumresource[vary,2,nreg] + median(novelwide[ , paste0(outcome, t,"allintermediate")])
    cumresource[vary,3,nreg] <- cumresource[vary,3,nreg] + median(novelwide[ , paste0(outcome, t, vary,"optimal")])
  }
}  


par(mar=c(6,6,6,1), mfcol=c(1,1), oma=c(0,0,0,0)) 
bres <- barplot(12*cbind(aperm(resource,c(1,3,2,4))["efficacy",c(2,1,3),,3], aperm(resource,c(1,3,2,4))["companion",c(2,1,3),,3],aperm(resource,c(1,3,2,4))["duration",c(2,1,3),,3]), beside = FALSE, 
                space=c(0.75,0.25,0.25), cex.lab=1.2, main="Treatment provided", ylim=c(0,500),
                col=rev(blues), ylab="Patient-months of treatment\nin year 10, by regimen", 
                names.arg=(c("minimal\n(lowest)","","optimal\n(highest)", 
                             "minimal\n(highest)","","optimal\n(lowest)",
                             "minimal\n(longest)","","optimal\n(shortest)")))
legend(x = bres[1], y=max(12*rowSums(resource["duration",,1:3,3]))+20, xjust=0.5, yjust=0, fill=blues,
       c("Novel DS regimen","MDR regimen","Standard DS regimen")[c(1,3,2)], xpd=NA)
text(bres[c(2,5,8)], -6*12, c("% Durably Cured", "Baseline novel-\nregimen resistance","Duration"), cex=1.1, pos=1, xpd=NA)


# diagnostic and other resources

tests <- array(0,dim=c( length(elementnames) , 3 , 3, 5 ));
dimnames(tests) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "test"=c("Novel regimen DSTs performed","Rifampin DSTs performed","Treatment courses initiated"),"q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in elementnames)  for (ntest in 1:3)
{ test <- c("nDSTs","rDSTs","dxs")[ntest]
  tests[vary,1,ntest,] <- quantile((1-novelwide$initialloss_s)*novelwide[ , paste0(test, "10", vary,"minimal")], c(0.025,0.25,0.5,0.75,0.975))
  tests[vary,2,ntest,] <- quantile((1-novelwide$initialloss_s)*novelwide[ , paste0(test, "10allintermediate")], c(0.025,0.25,0.5,0.75,0.975))
  tests[vary,3,ntest,] <- quantile((1-novelwide$initialloss_s)*novelwide[ , paste0(test, "10", vary,"optimal")], c(0.025,0.25,0.5,0.75,0.975))
}  

# par(mar=c(5,5.5,1,0),oma=c(0.5,0.5,0.5,0.5), mfrow=c(1,1))
# bup <- barplot(c(aperm(tests,c(1,3,2,4))["efficacy",3,,3], aperm(tests,c(1,3,2,4))["duration",3,,3], aperm(tests,c(1,3,2,4))["companion",3,,3]), beside = TRUE, 
#                space=c(0.75,0.25,0.25), ylim=c(0,100), cex.lab=1.2,
#                col=c("gray30","gray60","gray90"), ylab="Total treatment courses\ninitiated, year 10", xlab="")
# text(bup[c(2,5,8)], -20, elementlabels[2:4], cex=1, pos=1, xpd=NA)
# arrows(bup, aperm(tests, c(2,1,3,4))[,2:4,3,"0.25"], bup, aperm(tests, c(2,1,3,4))[,2:4,3,"0.75"], angle=90, code=3, length=0.05)

par(mar=c(5,6,2,1), mfcol=c(1,1), oma=c(0,0,0,0)) 
bup <- barplot(c(aperm(tests,c(1,3,2,4))["efficacy",1,,3], aperm(tests,c(1,3,2,4))["companion",1,,3],aperm(tests,c(1,3,2,4))["duration",1,,3]), beside = TRUE, 
               space=c(0.75,0.25,0.25), ylim=c(0,70), cex.lab=1.2, main="Diagnostic testing",
               col=cols, ylab="Novel regimen DSTs\nperformed, year 10", xlab="")
text(bup[c(2,5,8)], -10, c("% Durably Cured", "Baseline novel-\nregimen resistance","Duration"), cex=1.1, pos=1, xpd=NA)
arrows(bup, aperm(tests, c(2,1,3,4))[,2:4,1,"0.25"], bup, aperm(tests, c(2,1,3,4))[,2:4,1,"0.75"], angle=90, code=3, length=0.05)


# DR regimen:
novelwide <- allnovelwide[["DRDSTall_"]]; drout <- alldrout[1:nrow(novelwide),]

cumresource <- array(0,dim=c( length(elementnames) , 3 , 3)); resource <- array(0,dim=c( length(elementnames) , 3 , 3,5)); 
dimnames(cumresource) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "reg"=c("Standard DS regimen","Standard DR regimen","Novel regimen")); dimnames(resource) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "reg"=c("First-line","Second-line","Novel DR"), "q"=c(0.025,0.25,0.5,0.75,0.975))

for (vary in elementnames) for (nreg in 1:3) 
{ 
  outcome <- c("rxtime_s","rxtime_r","rxtime_n", "dxs")[nreg]
  resource[vary,1,nreg,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"minimal")],c(0.025,0.25,0.5,0.75,0.975))
  resource[vary,2,nreg,] <- quantile(novelwide[ , paste0(outcome, "10allintermediate")],c(0.025,0.25,0.5,0.75,0.975))
  resource[vary,3,nreg,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"optimal")],c(0.025,0.25,0.5,0.75,0.975))
  for (t in 1:10)
  {
    cumresource[vary,1,nreg] <- cumresource[vary,1,nreg] + median(novelwide[ , paste0(outcome, t, vary,"minimal")])
    cumresource[vary,2,nreg] <- cumresource[vary,2,nreg] + median(novelwide[ , paste0(outcome, t,"allintermediate")])
    cumresource[vary,3,nreg] <- cumresource[vary,3,nreg] + median(novelwide[ , paste0(outcome, t, vary,"optimal")])
  }
}  

# par(mar=c(5,6,6,1), mfcol=c(1,1), oma=c(0,0,0,0)) 
# bres <- barplot(12*cbind(aperm(resource,c(1,3,2,4))["efficacy",1:3,,3], aperm(resource,c(1,3,2,4))["companion",1:3,,3],aperm(resource,c(1,3,2,4))["duration",1:3,,3]), beside = FALSE, 
#                 space=c(0.75,0.25,0.25), cex.lab=1.2, main="Treatment provided", ylim=c(0,500),
#                 col=rev(blues), ylab="Patient-months of treatment\nin year 10, by regimen")
# legend(x = bres[1], y=max(12*rowSums(resource["duration",,1:3,3]))+10, xjust=0.5, yjust=0, fill=blues,
#        c("Novel DR regimen","Standard DR regimen","DS regimen"), xpd=NA)
# text(bres[c(2,5,8)], -8*12, c("% Durably Cured", "Baseline novel-\nregimen resistance","Duration"), cex=1.1, pos=1, xpd=NA)
# 
# par(mar=c(5,6,6,1), mfcol=c(1,1), oma=c(0,0,0,0)) 
# bres <- barplot(cbind(aperm(resource,c(1,3,2,4))["duration",2:3,,3]), beside = FALSE, 
#                 space=c(0.75,0.25,0.25), cex.lab=1.2, main="DR-specific treatment provided", ylim=c(0,3),
#                 col=rev(blues[1:2]), ylab="Patient-months of treatment\nin year 10, by regimen", names.arg=c("","",""))
# legend(x = bres[2], y=3, xjust=0, yjust=0, fill=blues[1:2],
#        c("Novel DR regimen","Standard DR regimen"), xpd=NA)
# text(bres, -0.2, c("20 month\nnovel regimen","9 month\nnovel regimen","6 month\nnovel regimen"), cex=1.1, pos=1, xpd=NA)

par(mar=c(7,6,4,1), mfcol=c(1,1), oma=c(0,0,0,0)) 
bres <- barplot(12*cbind(aperm(resource,c(1,3,2,4))["efficacy",2:3,,3], aperm(resource,c(1,3,2,4))["companion",2:3,,3],aperm(resource,c(1,3,2,4))["duration",2:3,,3]), beside = FALSE, 
                space=c(0.75,0.25,0.25), cex.lab=1.2, main="DR-specific treatment provided", ylim=c(0,30),
                col=rev(blues), ylab="Patient-months of treatment\nin year 10, by regimen")
legend(x = bres[3], y=max(12*rowSums(resource["duration",,2:3,3]))-4, xjust=0.5, yjust=0, fill=blues[2:3],
       c("Novel DR regimen","Standard DR regimen"), xpd=NA)
text(bres[7:9], -4, c("20mo","9mo","6mo"), cex=0.9, pos=1, xpd=NA)
text(bres[c(2,5,8)], -6, c("% Durably Cured", "Baseline novel-\nregimen resistance","Duration"), cex=1.1, pos=1, xpd=NA)


# par(mar=c(5,6,6,1), mfcol=c(1,1), oma=c(0,0,0,0)) 
# bres <- barplot(cbind(aperm(cumresource,c(1,3,2))["duration",2:3,]), beside = FALSE, 
#                 space=c(0.75,0.25,0.25), cex.lab=1.2, main="DR-specific treatment provided", ylim=c(0,50),
#                 col=rev(blues[1:2]), ylab="Cumulative patient-months of treatment, by regimen", names.arg=c("","",""))
# legend(x = bres[1], y=4, xjust=0.5, yjust=0, fill=blues[1:2],
#        c("Novel DR regimen","Standard DR regimen"), xpd=NA)
# text(bres, -0.2, c("20 month\nnovel regimen","9 month\nnovel regimen","6 month\nnovel regimen"), cex=1.1, pos=1, xpd=NA)

# diagnostic and other resources

tests <- array(0,dim=c( length(elementnames) , 3 , 3, 5 ));
dimnames(tests) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "test"=c("Novel regimen DSTs performed","Rifampin DSTs performed","Treatment courses initiated"),"q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in elementnames)  for (ntest in 1:3)
{ test <- c("nDSTs","rDSTs","dxs")[ntest]
  tests[vary,1,ntest,] <- quantile((1-novelwide$initialloss_s)*novelwide[ , paste0(test, "10", vary,"minimal")], c(0.025,0.25,0.5,0.75,0.975))
  tests[vary,2,ntest,] <- quantile((1-novelwide$initialloss_s)*novelwide[ , paste0(test, "10allintermediate")], c(0.025,0.25,0.5,0.75,0.975))
  tests[vary,3,ntest,] <- quantile((1-novelwide$initialloss_s)*novelwide[ , paste0(test, "10", vary,"optimal")], c(0.025,0.25,0.5,0.75,0.975))
}  

# par(mar=c(5,5.5,1,0),oma=c(0.5,0.5,0.5,0.5), mfrow=c(1,1))
# bup <- barplot(c(aperm(tests,c(1,3,2,4))["efficacy",3,,3], aperm(tests,c(1,3,2,4))["duration",3,,3], aperm(tests,c(1,3,2,4))["companion",3,,3]), beside = TRUE, 
#                space=c(0.75,0.25,0.25), ylim=c(0,100), cex.lab=1.2,
#                col=c("gray30","gray60","gray90"), ylab="Total treatment courses\ninitiated, year 10", xlab="")
# text(bup[c(2,5,8)], -20, elementlabels[2:4], cex=1, pos=1, xpd=NA)
# arrows(bup, aperm(tests, c(2,1,3,4))[,2:4,3,"0.25"], bup, aperm(tests, c(2,1,3,4))[,2:4,3,"0.75"], angle=90, code=3, length=0.05)

par(mar=c(5,6,2,1), mfcol=c(1,1), oma=c(0,0,0,0)) 
bup <- barplot(c(aperm(tests,c(1,3,2,4))["efficacy",1,,3], aperm(tests,c(1,3,2,4))["companion",1,,3],aperm(tests,c(1,3,2,4))["duration",1,,3]), beside = TRUE, 
               space=c(0.75,0.25,0.25), ylim=c(0,70), cex.lab=1.2, main="Diagnostic testing",
               col=cols, ylab="Novel regimen DSTs\nperformed, year 10", xlab="")
text(bup[c(2,5,8)], -15, c("% Durably Cured", "Baseline novel-regimen resistance","Duration"), cex=1.1, pos=1, xpd=NA)
arrows(bup, aperm(tests, c(2,1,3,4))[,2:4,1,"0.25"], bup, aperm(tests, c(2,1,3,4))[,2:4,1,"0.75"], angle=90, code=3, length=0.05)



#######
novelwide <- allnovelwide[["DSDSTall_rDSTall."]]; drout <- droutds <- alldrDST[1:nrow(novelwide),]

outcomes <- c("tbdeaths","inc","relapses") 

down <- array(0,dim=c( length(elementnames) , 3 , length(outcomes), 5 )); 
dimnames(down) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), outcome=outcomes, "q"=c(0.025,0.25,0.5,0.75,0.975))

for (vary in elementnames) for (nout in 1:length(outcomes))
{ 
  outcome <- outcomes[nout]
  down[vary,1,nout,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  down[vary,2,nout,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  down[vary,3,nout,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}  

bdown <- barplot(100*c(down["efficacy",,1,3], down["duration",,1,3], down["companion",,1,3]), beside = TRUE, 
                 space=c(0.75,0.25,0.25), ylim=c(-25,5), cex.lab=1.2,
                 col=cols, ylab="% reduction in TB deaths, year 10", xlab="", main="Mortality impact")
text(bup[c(2,5,8)], -28, elementlabels[2:4], cex=0.9, pos=1, xpd=NA)
# mtext("Varied TRP element (All others held fixed at intermediate level)", side=1, line=5, cex=0.7)
arrows(bdown, 100*c(t(down[2:4,,1,"0.25"])), bdown, 100*c(t(down[2:4,,1,"0.75"])), angle=90, code=3, length=0.05)

# tinc  <- array(0,dim=c( length(elementnames) , 3 , 5 )); 
# dimnames(tinc) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
# 
# for (vary in elementnames)
# { 
#   tinc[vary,1,] <- quantile((novelwide[ , paste0("inc", "10", vary,"minimal")] + novelwide[ , paste0("relapses", "10", vary,"minimal")] - 
#                                (drout[ , paste0("inc","10")]  + drout[ , paste0("relapses","10")] ))/
#                               (drout[ , paste0("inc","10")] + drout[ , paste0("relapses","10")] ), c(0.025,0.25,0.5,0.75,0.975))
#   tinc[vary,2,] <- quantile((novelwide[ , paste0("inc", "10allintermediate")] + novelwide[ , paste0("relapses", "10allintermediate")] - 
#                                (drout[ , paste0("inc","10")] + drout[ , paste0("relapses","10")] ))/
#                               (drout[ , paste0("inc","10")] + drout[ , paste0("relapses","10")]), c(0.025,0.25,0.5,0.75,0.975))
#   tinc[vary,3,] <- quantile((novelwide[ , paste0("inc", "10", vary,"optimal")] + novelwide[ , paste0("relapses", "10", vary,"optimal")] - 
#                                (drout[ , paste0("inc","10")]  + drout[ , paste0("relapses","10")] ))/
#                               (drout[ , paste0("inc","10")] + drout[ , paste0("relapses","10")] ), c(0.025,0.25,0.5,0.75,0.975))  
# }
# 
# btinc <- barplot(100*c(tinc["efficacy",,3], tinc["duration",,3], tinc["companion",,3]), beside = TRUE, 
#                  space=c(0.75,0.25,0.25), ylim=c(-15,2), cex.lab=1.2,
#                  col=cols, ylab="% reduction in TB incidence, year 10", xlab="")
# text(btinc[c(2,5,8)], -16.5, elementlabels[2:4], cex=1, pos=1, xpd=NA)
# arrows(btinc, 100*c(t(tinc[2:4,,"0.025"])), btinc, 100*c(t(tinc[2:4,,"0.975"])), angle=90, code=3, length=0.05)



# sensitivity analysis: ltfu at 2 months
ltfu2mo <- read.csv("TRPwideoutput_DSDSTall_ltfu2mo.rDSTall.India_20160201.1.csv", header=TRUE)
drltfu <- screendrout("DRcalibration_ltfu2mo.rDSTall.India_20160201.1.csv", tolerance=1.5)
outcome <- "tbdeaths"
down <- array(0,dim=c( length(elementnames) , 3 , 5 )); 
dimnames(down) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in elementnames)
{ down[vary,1,] <- quantile((ltfu2mo[ , paste0(outcome, "10", vary,"minimal")] - drltfu[ , paste0(outcome,"10")] )/
                              drltfu[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  down[vary,2,] <- quantile((ltfu2mo[ , paste0(outcome, "10allintermediate")] - drltfu[ , paste0(outcome,"10")] )/
                              drltfu[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  down[vary,3,] <- quantile((ltfu2mo[ , paste0(outcome, "10", vary,"optimal")] - drltfu[ , paste0(outcome,"10")] )/
                              drltfu[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}  

par(mar=c(2,4,3,1), mfrow=c(1,1), oma=c(0,0,0,0))
bpctdown <- barplot(height = 100*aperm(down, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="% reduction (median [IQR])",
                    main="With all losses to follow up occurring at 2 months\n(to maximize the impact of regimen duration)",
                    xlab="", las=2, cex.lab=1, 
                    ylim=c(150*final_pctdown[1,3,"0.5"], 2), yaxt='n',
#                     legend=c("minimal","intermediate","optimal"),
#                     args.legend=list(title="Level of varied element(s)", x=25,y=-15,cex=0.9),
                    space=c(0,0,0,1.5,0,0,rep(c(0.5,0,0),5)),
                    col=cols, names.arg=rep("", length(elementnames)))
axis(2, at=seq(-25,0,by=5), labels=paste0(seq(-25,0,by=5),"%"),las=2,cex.axis=0.8)
# text((bpctdown+0.4)[,2], -0.5, dslabels[4:6] ,cex=0.9, pos=2, srt=90, col="black", font=2) 
text(colMeans(bpctdown) ,-23, elementlabels, cex=0.9, pos=1, srt=90, font=1, xpd=NA)
mtext("Varied TRP element(s)", side=1, line=1, cex=0.9)
arrows(bpctdown, aperm(100*down, c(2,1,3))[,,"0.25"], bpctdown, aperm(100*down, c(2,1,3))[,,"0.75"], angle=90, code=3, length=0.05, xpd=NA)
# arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)


efracs <- (down[2:7,3,3]-down[2:7,1,3] )/sum(down[2:7,3,3]-down[2:7,1,3] )
xp <- numeric(6); for (i in 1:6) xp[i] <- sum(efracs[1:i])
xp2 <- c(0,xp[-6]) + (xp - c(0,xp[-6]))/2

ffracs <- (final_pctdown[2:7,3,3]-final_pctdown[2:7,1,3] )/sum(final_pctdown[2:7,3,3]-final_pctdown[2:7,1,3] )
fp <- numeric(6); for (i in 1:6) fp[i] <- sum(ffracs[1:i])
fp2 <- c(0,fp[-6]) + (fp - c(0,fp[-6]))/2

par(mfrow=c(1,1), mar=c(3,9,3,1), oma=c(0,0,0,0)) 
b <- barplot(array(c(efracs, ffracs), dim=c(6,2)), horiz = TRUE, beside=FALSE, las=2, 
             names.arg=c("With evenly-\ndistributed loss\nto follow up","With all losses to\nfollow up occurring\nat 2 months"), font=2,
             legend.text=elementnames[-1], col=rainbow(6), args.legend=list(x="left"))
mtext("Relative contributions of specified characteristics\nto novel regimen's impact", side=3, font=2, cex=1.1)
text(rbind(xp2,fp2)+0.005,b+0.1*(rep(c(0,0,1,1),3)),paste0(100*round(rbind(efracs,ffracs),2),"%"))


# sensitivity analysis: PRCCs for duration
library("sensitivity")
prcc <- list()
prcc$duration <- pcc(X = novelwide1[,40:73], y= (novelwide1[,"tbdeaths10durationminimal"] - novelwide1[,"tbdeaths10durationoptimal"] )/ (novelwide1[,"tbdeaths10allminimal"] - novelwide1[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$durationres <- pcc(X = novelwide3[,40:73], y= (novelwide3[,"tbdeaths10durationminimal"] - novelwide3[,"tbdeaths10durationoptimal"] )/ (novelwide3[,"tbdeaths10allminimal"] - novelwide3[,"tbdeaths10alloptimal"] ), rank=TRUE)
cbind(rownames(prcc$duration$PRCC), prcc$duration$PRCC)[rev(order(abs(prcc$duration$PRCC))),]

# plot PRCCs:
o <- rev(rev(order(abs(prcc$duration$PRCC$original)))[1:10])
ores <- rev(rev(order(abs(prcc$durationres$PRCC$original)))[1:10])
display <- c(ores[!(ores %in% o)],o)

par(mfrow=c(1,2),mar=c(4,2,2,1), oma=c(0,25,2,0))
b <- barplot(prcc$duration$PRCC$original[display], horiz=TRUE, beside=TRUE, space=c(0.5,0),
             names.arg=longparamnames[display], 
             axes=FALSE,axisnames=FALSE, xlim=c(-1,1),
             las=1, cex.names=0.8, cex.axis=0.8,cex.lab=0.8, main="Novel DS regimen", cex.main=1 )
axis(1, at=0.2*c(-4:4),cex.axis=0.8,cex.lab=0.9)
axis(2, labels=longparamnames[display],at=b, las=1, cex.axis=0.8, cex.lab=2, outer = TRUE, tick = FALSE, xpd=NA)
mtext("Parameter (varied +-50% around median estimate", side=2,line=24, outer=TRUE)
b <- barplot(prcc$durationres$PRCC$original[display], horiz=TRUE, beside=TRUE, space=c(0.5,0),
             names.arg=longparamnames[display], 
             axes=FALSE,axisnames=FALSE, xlim=c(-1,1),
             las=1, cex.names=0.8, cex.axis=0.8,cex.lab=0.8, main="Novel DR regimen", cex.main=1 )
axis(1, at=0.2*c(-4:4),cex.axis=0.8,cex.lab=0.9)
mtext("Partial rank correlation of model parameter with the mortality impact of regimen duration", side=3,cex=1.2,outer=TRUE)

# Note that in old (non rDSTall) DS scenario,
# Novel DS regimen results are highly sensitive to the fraction of DR TB being correctly diagnosed, 
# and to outcomes of DR TB on both first and second line therapy, because we've assumed the novel regimen works for DS TB too even it's not the intended target. 



###############
#Companion drug resistance

novelwide <- novelwide2 <- allnovelwide[["DSDSTnone_rDSTall."]]; drout <- alldrDST[1:nrow(novelwide),]

final_pctdown2 <- array(0,dim=c( length(elementnames) , 3 , 5 ));
dimnames(final_pctdown2) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))

outcome <- "tbdeaths"

for (vary in elementnames) 
{ 
  final_pctdown2[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  final_pctdown2[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  final_pctdown2[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}  

par(mar=c(5,5,3,1), mfrow=c(1,2), oma=c(0,0,0,0))

bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,c(4),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality (median [IQR])", 
                    xlab="Baseline prevalence of\ncompanion-drug resistance", names.arg=c("10%","3%","0%"), cex.lab=1, main="With DST for\nnovel DS regimen",
                    ylim=c(-15,2), col=cols)
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,4,"0.25"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,4,"0.75"], angle=90, code=3, length=0.05)


bpctdown <- barplot(height = 100*aperm(final_pctdown2, c(2,1,3))[,c(4),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality (median [IQR])", 
                    xlab="Baseline prevalence of\ncompanion-drug resistance", names.arg=c("10%","3%","0%"), cex.lab=1, main="Without DST for\nnovel DS regimen",
                    ylim=c(-15,2), col=cols)
arrows(bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,4,"0.25"], bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,4,"0.75"], angle=90, code=3, length=0.05)



# and same for DR regimen:# 
novelwide <- novelwide2 <- allnovelwide[["DRDSTnone_"]]; drout <- alldrout[1:nrow(novelwide),]

r_pctdown2 <- array(0,dim=c( length(elementnames) , 3 , 5 ));
dimnames(r_pctdown2) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))

outcome <- "rrdeaths"

for (vary in elementnames) 
{ 
  r_pctdown2[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  r_pctdown2[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  r_pctdown2[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}  

par(mar=c(5,5,3,1), mfrow=c(1,2), oma=c(0,0,0,0))

bpctdown <- barplot(height = 100*aperm(rr_pctdown, c(2,1,3))[,c(4),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 rifampin-resistant\nTB mortality (median [IQR])", 
                    xlab="Baseline prevalence of\ncompanion-drug resistance", names.arg=c("15%","5%","0%"), cex.lab=1, 
                    main="With DST for\nnovel DR regimen",
                    ylim=c(-40,2), col=cols)
arrows(bpctdown, aperm(100*rr_pctdown, c(2,1,3))[,4,"0.25"], bpctdown, aperm(100*rr_pctdown, c(2,1,3))[,4,"0.75"], angle=90, code=3, length=0.05)

bpctdown <- barplot(height = 100*aperm(r_pctdown2, c(2,1,3))[,c(4),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 rifampin-resistant\nTB mortality(median [IQR])", 
                    xlab="Baseline prevalence of\ncompanion-drug resistance", names.arg=c("15%","5%","0%"), cex.lab=1, 
                    main="Without DST for\nnovel DR regimen",
                    ylim=c(-40,2), col=cols)
arrows(bpctdown, aperm(100*r_pctdown2, c(2,1,3))[,4,"0.25"], bpctdown, aperm(100*r_pctdown2, c(2,1,3))[,4,"0.75"], angle=90, code=3, length=0.05)


#sensitivity
novelwide <- allnovelwide[["DSDSTall_rDSTall."]]; drout <- alldrDST[1:nrow(novelwide),]
prcc$companion <- pcc(X = novelwide[,40:73], y= (novelwide[,"tbdeaths10companionminimal"] - novelwide[,"tbdeaths10companionoptimal"] )/ (novelwide[,"tbdeaths10allminimal"] - novelwide[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$companionres <- pcc(X = novelwide[,40:73], y= ( novelwide[,"cprev10companionminimal"] / novelwide[,"cprev10companionoptimal"] ), rank=TRUE)
cbind(rownames(prcc$companion$PRCC), prcc$companion$PRCC)[rev(order(abs(prcc$companion$PRCC))),]

novelwide <- allnovelwide[["DRDSTall_"]]; drout <- alldrout[1:nrow(novelwide),]
prcc$companion_dr <- pcc(X = novelwide[,40:73], y= (novelwide[,"tbdeaths10companionminimal"] - novelwide[,"tbdeaths10companionoptimal"] )/ (novelwide[,"tbdeaths10allminimal"] - novelwide[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$companionres_dr <- pcc(X = novelwide[,40:73], y= ( novelwide[,"novelprev10companionminimal"] / novelwide[,"novelprev10companionoptimal"] ), rank=TRUE)
# summary: correlated with parameters that increase the number of treatments (and to a lesser extent, that reduce the new dr transmissions to dilute the remaining c cases) per incident case

novelwide <- allnovelwide[["DSDSTnone_rDSTall."]]; drout <- alldrDST[1:nrow(novelwide),]
prcc$companionnodst <- pcc(X = novelwide[,40:73], y= (novelwide[,"tbdeaths10companionminimal"] - novelwide[,"tbdeaths10companionoptimal"] )/ (novelwide[,"tbdeaths10allminimal"] - novelwide[,"tbdeaths10alloptimal"] ), rank=TRUE)

novelwide <- allnovelwide[["DRDSTnone_"]]; drout <- alldrout[1:nrow(novelwide),]
prcc$companionnodst_dr <- pcc(X = novelwide[,40:73], y= (novelwide[,"tbdeaths10companionminimal"] - novelwide[,"tbdeaths10companionoptimal"] )/ (novelwide[,"tbdeaths10allminimal"] - novelwide[,"tbdeaths10alloptimal"] ), rank=TRUE)

# Plot PRCCs with DST:
o <- rev(rev(order(abs(prcc$companion$PRCC$original)))[1:10])
ores <- rev(rev(order(abs(prcc$companion_dr$PRCC$original)))[1:10])
display <- c(ores[!(ores %in% o)],o)

par(mfrow=c(1,2),mar=c(4,0,2,2), oma=c(0,25,3,0))
b <- barplot(prcc$companion$PRCC$original[display], horiz=TRUE, beside=TRUE, space=c(0.5,0),
             names.arg=longparamnames[display], 
             axes=FALSE,axisnames=FALSE, xlim=c(-1,1),
             las=1, cex.names=0.8, cex.axis=0.8,cex.lab=0.8, main="Novel DS regimen", cex.main=1 )
axis(1, at=0.2*c(-4:4),cex.axis=0.8,cex.lab=0.9)
axis(2, labels=longparamnames[display],at=b, las=1, cex.axis=0.8, cex.lab=2, outer = TRUE, tick = FALSE, xpd=NA)
mtext("Parameter (varied +-50% around median estimate)", side=2, line=24, outer=TRUE)
b <- barplot(prcc$companion_dr$PRCC$original[display], horiz=TRUE, beside=TRUE, space=c(0.5,0),
             names.arg=longparamnames[display], 
             axes=FALSE,axisnames=FALSE, xlim=c(-1,1),
             las=1, cex.names=0.8, cex.axis=0.8,cex.lab=0.8, main="Novel DR regimen", cex.main=1 )
axis(1, at=0.2*c(-4:4),cex.axis=0.8,cex.lab=0.9)
mtext("Partial rank correlation of model parameter with the\nmortality impact of companion-drug resistance prevalence, with DST", side=3,cex=1.2,outer=TRUE)


# plot PRCCs for no DST:
o <- rev(rev(order(abs(prcc$companionnodst$PRCC$original)))[1:10])
ores <- rev(rev(order(abs(prcc$companionnodst_dr$PRCC$original)))[1:10])
display <- c(ores[!(ores %in% o)],o)

par(mfrow=c(1,2),mar=c(4,0,2,2), oma=c(0,25,3,0))
b <- barplot(prcc$companionnodst$PRCC$original[display], horiz=TRUE, beside=TRUE, space=c(0.5,0),
             names.arg=longparamnames[display], 
             axes=FALSE,axisnames=FALSE, xlim=c(-1,1),
             las=1, cex.names=0.8, cex.axis=0.8,cex.lab=0.8, main="Novel DS regimen", cex.main=1 )
axis(1, at=0.2*c(-4:4),cex.axis=0.8,cex.lab=0.9)
axis(2, labels=longparamnames[display],at=b, las=1, cex.axis=0.8, cex.lab=2, outer = TRUE, tick = FALSE, xpd=NA)
mtext("Parameter (varied +-50% around median estimate)", side=2, line=24, outer=TRUE)
b <- barplot(prcc$companionnodst_dr$PRCC$original[display], horiz=TRUE, beside=TRUE, space=c(0.5,0),
             names.arg=longparamnames[display], 
             axes=FALSE,axisnames=FALSE, xlim=c(-1,1),
             las=1, cex.names=0.8, cex.axis=0.8,cex.lab=0.8, main="Novel DR regimen", cex.main=1 )
axis(1, at=0.2*c(-4:4),cex.axis=0.8,cex.lab=0.9)
mtext("Partial rank correlation of model parameter with the\nmortality impact of companion-drug resistance prevalence, without DST", side=3,cex=1.2,outer=TRUE)


###########################
## Barrier to resistance ##

#DS:
par(mar=c(5,5,3,1), mfrow=c(1,2), oma=c(0,0,0,0))

bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,c(5),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality (median [IQR])", 
                    xlab="Per-treatment probability\nof acquired novel resistance", names.arg=c("1/20","1/125","0"), cex.lab=1, main="With DST for\nnovel DS regimen",
                    ylim=c(-15,2), col=cols)
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,5,"0.25"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,5,"0.75"], angle=90, code=3, length=0.05)


bpctdown <- barplot(height = 100*aperm(final_pctdown2, c(2,1,3))[,c(5),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality (median [IQR])", 
                    xlab="Per-treatment probability\nof acquired novel resistance", names.arg=c("1/20","1/125","0"), cex.lab=1, main="Without DST for\nnovel DS regimen",
                    ylim=c(-15,2), col=cols)
arrows(bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,5,"0.25"], bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,5,"0.75"], angle=90, code=3, length=0.05)


#DR:
par(mar=c(5,5,3,1), mfrow=c(1,2), oma=c(0,0,0,0))

bpctdown <- barplot(height = 100*aperm(rr_pctdown, c(2,1,3))[,c(5),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 rifampin-resistant\nTB mortality (median [IQR])", 
                    xlab="Per-treatment probability\nof acquired novel resistance", names.arg=c("1/10","1/20","1/125"), cex.lab=1, 
                    main="With DST for\nnovel DR regimen",
                    ylim=c(-40,2), col=cols)
arrows(bpctdown, aperm(100*rr_pctdown, c(2,1,3))[,5,"0.25"], bpctdown, aperm(100*rr_pctdown, c(2,1,3))[,5,"0.75"], angle=90, code=3, length=0.05)

bpctdown <- barplot(height = 100*aperm(r_pctdown2, c(2,1,3))[,c(5),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 rifampin-resistant\nTB mortality(median [IQR])", 
                    xlab="Per-treatment probability\nof acquired novel resistance", names.arg=c("1/10","1/20","1/125"), cex.lab=1, 
                    main="Without DST for\nnovel DR regimen",
                    ylim=c(-40,2), col=cols)
arrows(bpctdown, aperm(100*r_pctdown2, c(2,1,3))[,5,"0.25"], bpctdown, aperm(100*r_pctdown2, c(2,1,3))[,5,"0.75"], angle=90, code=3, length=0.05)


#sensitivity: PRCCs for barrier to resistance
novelwide <- allnovelwide[["DSDSTall_rDSTall."]]; drout <- alldrDST[1:nrow(novelwide),]
prcc$barrier <- pcc(X = novelwide[,40:73], y= (novelwide[,"tbdeaths10barrierminimal"] - novelwide[,"tbdeaths10barrieroptimal"] )/ (novelwide[,"tbdeaths10allminimal"] - novelwide[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$barrierres <- pcc(X = novelwide[,40:73], y= ( novelwide[,"cprev10barrierminimal"] / novelwide[,"cprev10barrieroptimal"] ), rank=TRUE)
cbind(rownames(prcc$barrier$PRCC), prcc$barrier$PRCC)[rev(order(abs(prcc$barrier$PRCC))),]

novelwide <- allnovelwide[["DRDSTall_"]]; drout <- alldrout[1:nrow(novelwide),]
prcc$barrier_dr <- pcc(X = novelwide[,40:73], y= (novelwide[,"tbdeaths10barrierminimal"] - novelwide[,"tbdeaths10barrieroptimal"] )/ (novelwide[,"tbdeaths10allminimal"] - novelwide[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$barrierres_dr <- pcc(X = novelwide[,40:73], y= ( novelwide[,"novelprev10barrierminimal"] / novelwide[,"novelprev10barrieroptimal"] ), rank=TRUE)
# summary: correlated with parameters that increase the number of treatments (and to a lesser extent, that reduce the new dr transmissions to dilute the remaining c cases) per incident case

novelwide <- allnovelwide[["DSDSTnone_rDSTall."]]; drout <- alldrDST[1:nrow(novelwide),]
prcc$barriernodst <- pcc(X = novelwide[,40:73], y= (novelwide[,"tbdeaths10barrierminimal"] - novelwide[,"tbdeaths10barrieroptimal"] )/ (novelwide[,"tbdeaths10allminimal"] - novelwide[,"tbdeaths10alloptimal"] ), rank=TRUE)

novelwide <- allnovelwide[["DRDSTnone_"]]; drout <- alldrout[1:nrow(novelwide),]
prcc$barriernodst_dr <- pcc(X = novelwide[,40:73], y= (novelwide[,"tbdeaths10barrierminimal"] - novelwide[,"tbdeaths10barrieroptimal"] )/ (novelwide[,"tbdeaths10allminimal"] - novelwide[,"tbdeaths10alloptimal"] ), rank=TRUE)

# Plot PRCCs with DST:
o <- rev(rev(order(abs(prcc$barrier$PRCC$original)))[1:10])
ores <- rev(rev(order(abs(prcc$barrier_dr$PRCC$original)))[1:10])
display <- c(ores[!(ores %in% o)],o)

par(mfrow=c(1,2),mar=c(4,0,2,2), oma=c(0,25,3,0))
b <- barplot(prcc$barrier$PRCC$original[display], horiz=TRUE, beside=TRUE, space=c(0.5,0),
             names.arg=longparamnames[display], 
             axes=FALSE,axisnames=FALSE, xlim=c(-1,1),
             las=1, cex.names=0.8, cex.axis=0.8,cex.lab=0.8, main="Novel DS regimen", cex.main=1 )
axis(1, at=0.2*c(-4:4),cex.axis=0.8,cex.lab=0.9)
axis(2, labels=longparamnames[display],at=b, las=1, cex.axis=0.8, cex.lab=2, outer = TRUE, tick = FALSE, xpd=NA)
mtext("Parameter (varied +-50% around median estimate)", side=2, line=24, outer=TRUE)
b <- barplot(prcc$barrier_dr$PRCC$original[display], horiz=TRUE, beside=TRUE, space=c(0.5,0),
             names.arg=longparamnames[display], 
             axes=FALSE,axisnames=FALSE, xlim=c(-1,1),
             las=1, cex.names=0.8, cex.axis=0.8,cex.lab=0.8, main="Novel DR regimen", cex.main=1 )
axis(1, at=0.2*c(-4:4),cex.axis=0.8,cex.lab=0.9)
mtext("Partial rank correlation of model parameter with the\nmortality impact of barrier to acquired novel-drug resistance, with DST", side=3,cex=1.2,outer=TRUE)


# plot PRCCs for no DST:
o <- rev(rev(order(abs(prcc$barriernodst$PRCC$original)))[1:10])
ores <- rev(rev(order(abs(prcc$barriernodst_dr$PRCC$original)))[1:10])
display <- c(ores[!(ores %in% o)],o)

par(mfrow=c(1,2),mar=c(4,0,2,2), oma=c(0,30,3,0))
b <- barplot(prcc$barriernodst$PRCC$original[display], horiz=TRUE, beside=TRUE, space=c(0.5,0),
             names.arg=longparamnames[display], 
             axes=FALSE,axisnames=FALSE, xlim=c(-1,1),
             las=1, cex.names=0.8, cex.axis=0.8,cex.lab=0.8, main="Novel DS regimen", cex.main=1 )
axis(1, at=0.2*c(-4:4),cex.axis=0.8,cex.lab=0.9)
axis(2, labels=longparamnames[display],at=b, las=1, cex.axis=0.8, cex.lab=2, outer = TRUE, tick = FALSE, xpd=NA)
mtext("Parameter (varied +-50% around median estimate)", side=2, line=29, outer=TRUE)
b <- barplot(prcc$barriernodst_dr$PRCC$original[display], horiz=TRUE, beside=TRUE, space=c(0.5,0),
             names.arg=longparamnames[display], 
             axes=FALSE,axisnames=FALSE, xlim=c(-1,1),
             las=1, cex.names=0.8, cex.axis=0.8,cex.lab=0.8, main="Novel DR regimen", cex.main=1 )
axis(1, at=0.2*c(-4:4),cex.axis=0.8,cex.lab=0.9)
mtext("Partial rank correlation of model parameter with the\nmortality impact of barrier to acquired novel-drug resistance, without DST", side=3,cex=1.2,outer=TRUE)


# Trends in resistance
## problem with main results: measuring total and drug-resistant prevalences differently; so ran separate "Resistance" simulations doing it correctly.

rtallynames <- readRDS("rtallynames_20160111.RDS")

RtrajDST <- read.csv("Resistance_DRDSTall_India_20160201.idr1.csv")
for (i in 2:5) RtrajDST <- rbind(RtrajDST, read.csv(paste0("Resistance_DRDSTall_India_20160201.idr",i,".csv")))
drout <- rbind(alldrDST[alldrDST$idr == 1,], alldrDST[alldrDST$idr == 2,], alldrDST[alldrDST$idr == 3,],
               alldrDST[alldrDST$idr == 4,], alldrDST[alldrDST$idr == 2,])

rtraj <- array(0,dim=c( 3, 3 , 5, 5, 11));
dimnames(rtraj) <- list("barrier"=levels, "companion"=levels, 
                         "drugs"=c("nprev","cprev","novelprev","rprev","rnovelprev"),"q"=c(0.025,0.25,0.5,0.75,0.975), "time"=0:10)
longdrugnames <- c("Novel drug","Companion drug","Novel and/or companion drug", "Rifampin", "Rifampin and (novel and/or companion) drug")

for (time in as.character(0:10)) for (barrier in levels) for (companion in levels)  for (ndrug in 1:5) 
{ 
  drug <- dimnames(rtraj)$drugs[ndrug]
  rtraj[barrier,companion,drug,,time] <- quantile(RtrajDST[ , paste0(drug, time, "companion",companion,"barrier",barrier)]/
                                                    RtrajDST[ , paste0("aprev", time, "companion",companion,"barrier",barrier)], c(0.025,0.25,0.5,0.75,0.975))
}  

par(mar=c(3,3,1,1),oma=c(0,4,5,0), mfrow=c(2,2))
for (b in c(3,1)) for (c in c(3,1)) 
{plot(0:10, rtraj[b,c,"novelprev",3,], type='l', lwd=2, ylim=c(0,0.2), xlab="",ylab="Prevalence of resistance", cex.lab=0.8, yaxt="n");
 if (b==3 & c==3) legend("topright",c("Any novel-regimen resistance","Companion-drug resistance", "Novel-drug resistance"), lwd=2,lty=1,col=c("black","green","purple"), cex=0.9)
 if (b==3) mtext(paste0(c("Unfavorable (10%)","Intermediate (3%)","No")[c]," baseline\ncompanion-drug resistance"),side=3, line =1, cex=0.8,font=2)
 if (c==3) mtext(paste0(c("Minimal","Intermediate","Optimal")[b]," barrier to resistance\n",c("(1/20 treatments","(1/125 treatments", "(no acquired resistance")[b],"\nif iniitally regimen- susceptible)"),side=2, cex=0.8,line=4, font=2)
#   mtext(paste0(c("Unfavorable (10%)","Optimal (3%)","No")[c]," baseline companion-resistance,\nand ",levels[b]," barrier to resistance\n(",c("1/20 treatments)","1/125 treatments)", "no acquired resistance if iniitally fully susceptible")[b]),
#         side=3,line=1,xpd=NA, cex=0.7)
  mtext("Year",side=1,line=2,cex=0.7)
  axis(side=2,at = c(seq(0,0.3,by=0.05)), labels = paste0(seq(0,30,by=5),"%"),las=2)
points(0:10, rtraj[b,c,"nprev",3,], type='l', lwd=2,col="purple")
points(0:10, rtraj[b,c,"cprev",3,], type='l', lwd=2,col="green")
}
mtext("With DST for novel regimen",side=3,outer=TRUE,line=3,cex=1.3)

Rtraj <- read.csv("Resistance_DRDSTnone_India_20160201.idr1.csv")
for (i in 2:5) Rtraj <- rbind(Rtraj, read.csv(paste0("Resistance_DRDSTnone_India_20160201.idr",i,".csv")))
drout <- rbind(alldrDST[alldrDST$idr == 1,], alldrDST[alldrDST$idr == 2,], alldrDST[alldrDST$idr == 3,],
               alldrDST[alldrDST$idr == 4,], alldrDST[alldrDST$idr == 2,])

rtraj <- array(0,dim=c( 3, 3 , 5, 5, 11));
dimnames(rtraj) <- list("barrier"=levels, "companion"=levels, 
                        "drugs"=c("nprev","cprev","novelprev","rprev","rnovelprev"),"q"=c(0.025,0.25,0.5,0.75,0.975), "time"=0:10)
longdrugnames <- c("Novel drug","Companion drug","Novel and/or companion drug", "Rifampin", "Rifampin and (novel and/or companion) drug")

for (time in as.character(0:10)) for (barrier in levels) for (companion in levels)  for (ndrug in 1:5) 
{ 
  drug <- dimnames(rtraj)$drugs[ndrug]
  rtraj[barrier,companion,drug,,time] <- quantile(Rtraj[ , paste0(drug, time, "companion",companion,"barrier",barrier)]/
                                                    Rtraj[ , paste0("aprev", time, "companion",companion,"barrier",barrier)], c(0.025,0.25,0.5,0.75,0.975))
}  

par(mar=c(3,3,1,1),oma=c(0,4,5,0), mfrow=c(2,2))
for (b in c(3,1)) for (c in c(3,1)) 
{plot(0:10, rtraj[b,c,"novelprev",3,], type='l', lwd=2, ylim=c(0,0.2), xlab="",ylab="Prevalence of resistance", cex.lab=0.8, yaxt="n");
  if (b==3 & c==3) legend("topright",c("Any novel-regimen resistance","Companion-drug resistance", "Novel-drug resistance"), lwd=2,lty=1,col=c("black","green","purple"), cex=0.9)
 if (b==3) mtext(paste0(c("Unfavorable (10%)","Intermediate (3%)","No")[c]," baseline\ncompanion-drug resistance"),side=3, line =1, cex=0.8,font=2)
 if (c==3) mtext(paste0(c("Minimal","Intermediate","Optimal")[b]," barrier to resistance\n",c("(1 in 20 acquire resistance","(1 in 125 acquire resistance", "(no acquired resistance")[b],"\nif iniitally regimen- susceptible)"),side=2, cex=0.8,line=4, font=2) #   mtext(paste0(c("Unfavorable (10%)","Optimal (3%)","No")[c]," baseline companion-resistance,\nand ",levels[b]," barrier to resistance\n(",c("1/20 treatments)","1/125 treatments)", "no acquired resistance if iniitally fully susceptible")[b]),
 #         side=3,line=1,xpd=NA, cex=0.7)
 mtext("Year",side=1,line=2,cex=0.7)
 axis(side=2,at = c(seq(0,0.3,by=0.05)), labels = paste0(seq(0,30,by=5),"%"),las=2)
 points(0:10, rtraj[b,c,"nprev",3,], type='l', lwd=2,col="purple")
 points(0:10, rtraj[b,c,"cprev",3,], type='l', lwd=2,col="green")
}
mtext("Without DST for novel regimen",side=3,outer=TRUE,line=3,cex=1.3)



# Projections without novel regimen #

par(mar=c(4,4,4,1))
rrfracs <- apply(alldrout[,paste0("rrfrac",c("","10"))], 2, quantile, c(0.025,0.25,0.5,0.75,0.975))
colnames(rrfracs) <- c("year 0", "year 10 without novel regimen")

rrfracsDST <- apply(alldrDST[,paste0("rrfrac",c("","10"))], 2, quantile, c(0.025,0.25,0.5,0.75,0.975))
colnames(rrfracsDST) <- c("year 0", "year 10 without novel regimen")


############################
# exclusions

# India exclusions

iexc <- rbind(read.csv("Exclusions_DRDSTall_India_20160111.idr1.csv",header=FALSE), read.csv("Exclusions_DRDSTall_India_20160111.idr2.csv",header=FALSE))
# colnames(iexc) <- colnames(read.csv("Exclusions_DSDSTall_rDSTall.India_20160111.idr1.csv"))
drout <- rbind(alldrout[alldrout$idr == 1,], alldrout[alldrout$idr == 2,])
outcome <- "rrdeaths"
idown <- array(0,dim=c(  3, 3, 2, 5 ));
dimnames(idown) <- list( "efflevel"=levels, "exclevel"=levels, "HIV"=c("HIV","nonHIV"), "q"=c(0.025,0.25,0.5,0.75,0.975))

for (efficacy in levels) for (exclusions in levels) for (H in c("HIV","nonHIV"))
{  idown[efficacy,exclusions, H,] <- quantile((iexc[ , paste0(outcome, "10", H,"exclusions",exclusions,"efficacy",efficacy)] - drout[ , paste0(outcome,"10")] )/
                               drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}

par(mar=c(5,5,1,1), mfrow=c(1,2), oma=c(0,0,2,0))

inoH <- barplot(height = 100*idown["intermediate",,"nonHIV","0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality (median [IQR])", 
                    xlab="HIV-unrelated exclusions\nfrom novel regimen", names.arg=c("11%","5%","0%"), cex.lab=1, main="",
                    ylim=c(-40,2), col=cols)
mtext("Novel DR regimen; intermediate efficacy; 5% HIV coprevalence", outer=TRUE, side=3, cex=1.4)
arrows(inoH, 100*idown["intermediate",,"nonHIV","0.25"], inoH, 100*idown["intermediate",,"nonHIV","0.75"], angle=90, code=3, length=0.05)


iH <- barplot(height = 100*idown["intermediate",,"HIV","0.5"], beside = TRUE, 
                  ylab="% reduction in year 10 TB mortality (median [IQR])", 
                  xlab="HIV-related exclusions\nfrom novel regimen", names.arg=c("100%","5%","0%"), cex.lab=1, main="",
                  ylim=c(-40,2), col=cols)
arrows(inoH, 100*idown["intermediate",,"HIV","0.25"], inoH, 100*idown["intermediate",,"HIV","0.75"], angle=90, code=3, length=0.05)


# South Africa exclusions

sexc <- rbind(read.csv("Exclusions_DRDSTall_SouthAfrica_20160111.idr1.csv"), read.csv("Exclusions_DRDSTall_SouthAfrica_20160111.idr2.csv"))
sdrout <- read.csv("DRcalibration_SouthAfrica_20160111.1.csv"); 
  sdrout <- sdrout[sdrout[,"rrinc"]/sdrout[,"inc"] > 1/tolerance*sdrout[,"targetdr"] & sdrout[,"rrinc"]/sdrout[,"inc"] < tolerance*sdrout[,"targetdr"], ];
  drout <- rbind(sdrout[sdrout$idr == 1,], sdrout[sdrout$idr == 2,])
outcome <- "rrdeaths"
sdown <- array(0,dim=c( 3, 3, 2, 5 ));
dimnames(sdown) <- list("efflevel"=levels, "exclevel"=levels, "HIV"=c("HIV","nonHIV"), "q"=c(0.025,0.25,0.5,0.75,0.975))

for (efficacy in levels) for (exclusions in levels) for (H in c("HIV","nonHIV"))
{  sdown[efficacy,exclusions, H,] <- quantile((sexc[ , paste0(outcome, "10", H,"exclusions",exclusions,"efficacy",efficacy)] - drout[ , paste0(outcome,"10")] )/
                                                     drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}

par(mar=c(5,5,1,1), mfrow=c(1,2), oma=c(0,0,2,0))

snoH <- barplot(height = 100*sdown["intermediate",,"nonHIV","0.5"], beside = TRUE, 
                ylab="% reduction in year 10 TB mortality (median [IQR])", 
                xlab="HIV-unrelated exclusions\nfrom novel regimen", names.arg=c("11%","5%","0%"), cex.lab=1, main="",
                ylim=c(-60,2), col=cols)
mtext("Novel DR regimen; intermediate efficacy; 60% HIV coprevalence", outer=TRUE, side=3, cex=1.4)
arrows(snoH, 100*sdown["intermediate",,"nonHIV","0.25"], snoH, 100*sdown["intermediate",,"nonHIV","0.75"], angle=90, code=3, length=0.05)


sH <- barplot(height = 100*sdown["intermediate",,"HIV","0.5"], beside = TRUE, 
              ylab="% reduction in year 10 TB mortality (median [IQR])", 
              xlab="HIV-related exclusions\nfrom novel regimen", names.arg=c("100%","5%","0%"), cex.lab=1, main="",
              ylim=c(-60,2), col=cols)
arrows(sH, 100*sdown["intermediate",,"HIV","0.25"], sH, 100*sdown["intermediate",,"HIV","0.75"], angle=90, code=3, length=0.05)

# # synergy with maximal efficacy
# 
# par(mfrow=c(1,1), mar=c(9,3,4,1))
# bpctdown <- barplot(height = 100*aperm(edown, c(2,1,3))[,,"0.5"], beside = TRUE, 
#                     col=c("gray30","gray60","gray90"), ylim=100*min(edown[,,"0.5"])*c(1.5,-0.1),
#                     names.arg=elementlabels[2:7], cex.names=0.7,
#                     main="Exploring synergy with efficacy:\nNovel regimen impact with maximal % durably cured throughout\n(but for other characteristics still intermediate level)")
# arrows(bpctdown, aperm(100*edown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*edown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)
# 
# text(mean(bpctdown), -32.5, "Fraction of total variation in novel regimen's impact attributable to variation in characteristic:
#      \nWith novel regimen efficacy at maximal level (99% cure) throughout", pos=1, xpd=NA)
# text(bpctdown[2,], -35.5, paste0(round((edown[,1,3]-edown[,3,3] )/(median(novelwide[,"tbdeaths10allminimal"]/drout[,"tbdeaths10"])-median(novelwide[,"tbdeaths10alloptimal"]/drout[,"tbdeaths10"]))*100, digits=1), "%") , pos=1, xpd=NA)
# text(mean(bpctdown), -37, "With novel regimen efficacy at intermediate level (97% cure) throughout", pos=1, xpd=NA)
# text(bpctdown[2,], -38, paste0(round(((final_pctdown[,3,3]-final_pctdown[,1,3] )/(final_pctdown[1,3,3]-final_pctdown[1,1,3] ))[2:7]*100, digits=1), "%") , pos=1, xpd=NA)

# synergy: impact of exclusion variation when efficacy intermediate vs when efficacy optimal: 

inoHlow <- barplot(height = 100*idown["minimal",,"nonHIV","0.5"], beside = TRUE, 
                ylab="% reduction in year 10 TB mortality (median [IQR])", 
                xlab="HIV-unrelated exclusions\nfrom novel regimen", names.arg=c("11%","5%","0%"), cex.lab=1, main="",
                ylim=c(-15,2), col=cols)
mtext("Minimal efficacy", outer=FALSE, side=3, cex=1.2)
arrows(inoHlow, 100*idown["minimal",,"nonHIV","0.25"], inoHlow, 100*idown["minimal",,"nonHIV","0.75"], angle=90, code=3, length=0.05)

# iHlow <- barplot(height = 100*idown["minimal",,"HIV","0.5"], beside = TRUE, 
#               ylab="% reduction in year 10 TB mortality (median [IQR])", 
#               xlab="HIV-related exclusions\nfrom novel regimen", names.arg=c("100%","5%","0%"), cex.lab=1, main="",
#               ylim=c(-15,2), col=cols)
# arrows(iHlow, 100*idown["minimal",,"HIV","0.25"], iHlow, 100*idown["minimal",,"HIV","0.75"], angle=90, code=3, length=0.05)


inoHhigh <- barplot(height = 100*idown["optimal",,"nonHIV","0.5"], beside = TRUE, 
                   ylab="% reduction in year 10 TB mortality (median [IQR])", 
                   xlab="HIV-unrelated exclusions\nfrom novel regimen", names.arg=c("11%","5%","0%"), cex.lab=1, main="",
                   ylim=c(-30,2), col=cols)
mtext("Optimal efficacy", outer=FALSE, side=3, cex=1.2)
arrows(inoHhigh, 100*idown["optimal",,"nonHIV","0.25"], inoHhigh, 100*idown["optimal",,"nonHIV","0.75"], angle=90, code=3, length=0.05)

# 
# iHhigh <- barplot(height = 100*idown["optimal",,"HIV","0.5"], beside = TRUE, 
#                  ylab="% reduction in year 10 TB mortality (median [IQR])", 
#                  xlab="HIV-related exclusions\nfrom novel regimen", names.arg=c("100%","5%","0%"), cex.lab=1, main="",
#                  ylim=c(-30,2), col=cols)
# arrows(iHhigh, 100*idown["optimal",,"HIV","0.25"], iHlow, 100*idown["optimal",,"HIV","0.75"], angle=90, code=3, length=0.05)



#####################
# Adherence #

par(mar=c(5,5,3,1), mfrow=c(1,1), oma=c(0,0,0,0))

bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,c(7),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality (median [IQR])", 
                    xlab="Improvement in adherence", names.arg=c("same as\nSOC","1.5% more\nfinish", "3% more\nfinish"), cex.lab=1,
                    ylim=c(-15,2), col=cols, main="Impact of\nimproved adherence")
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,7,"0.25"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,7,"0.75"], angle=90, code=3, length=0.05)


bpctdown <- barplot(height = 100*aperm(final_pctdown2, c(2,1,3))[,c(4),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality (median [IQR])", 
                    xlab="Baseline prevalence of\ncompanion-drug resistance", names.arg=c("10%","3%","0%"), cex.lab=1, main="Without DST for\nnovel DS regimen",
                    ylim=c(-15,2), col=cols)
arrows(bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,4,"0.25"], bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,4,"0.75"], angle=90, code=3, length=0.05)


par(mar=c(5,6,3,2), mfcol=c(1,1), oma=c(0,0,0,0)) 
adh <- barplot(12*cbind(aperm(resource,c(1,3,2,4))["efficacy",1:3,,3], aperm(resource,c(1,3,2,4))["duration",1:3,,3], aperm(resource,c(1,3,2,4))["companion",1:3,,3]), beside = FALSE, 
                space=c(0.75,0.25,0.25), cex.lab=1.2, main="Treatment costs",
                col=rev(blues), ylab="Patient-months of treatment\nin year 10, by regimen")
legend(x = bres[6], y=max(12*rowSums(resource["efficacy",,1:3,3]))+10, xjust=0, yjust=0, fill=blues,
       c("Novel DS regimen","Second-line regimen","First-line regimen"), xpd=NA)
text(bres[c(2,5,8)], -7*12, elementlabels[2:4], cex=0.9, pos=1, xpd=NA)
