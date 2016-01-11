source("TPPmat.R")

currenttag <- "India_20160105"; tolerance <- 1.5; location=""

source("displayfunction.R")


###############################
# Results overview #

novelwide <- allnovelwide[["DSDSTall"]]; drout <- alldrout[1:nrow(novelwide),]

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
plot(0:10, traj[,"minimal","0.5"], ylim=c(0,targetepis[[targetepi]][1]/4), type='l', col="red", xlab="Years after introduction", ylab="Annual TB mortality")
points(0:10, traj[,"intermediate","0.5"], type='l', col='orange')
points(0:10, traj[,"optimal","0.5"], type='l', col='green')
points(c(0,10), c(median(drout[,"tbdeaths"]),median(drout[,"tbdeaths10"])), col='black', type='l') 
legend("bottomleft", legend= c("Continued current standard", "Minimal novel DS regimen", "Intermediate novel DS regimen", "Optimal novel DS regimen"),
       col=c("black","red","orange","green"), lty=1)

par(mar=c(3,5,4,3), mfrow=c(1,1))
bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="Percent reduction in year 10 TB mortality with novel regimen,\n compared to projection under current standard of care", 
                    xlab="", las=2, cex.lab=1, 
                    ylim=c(150*final_pctdown[1,3,"0.5"], 8),
                    legend=c("minimal","intermediate","optimal"), args.legend=list(title="Level of varied element(s)",x=bpctdown[2,7], y=110*final_pctdown[1,3,"0.5"],cex=0.8),
                    main="", #paste0("Novel regimen for DS TB, with universal DST"), cex.main=1,
                    col=cols, names.arg=rep("", length(elementnames)))
text(bpctdown+0.4, -0.5, dslabels ,cex=0.9, pos=2, srt=90, col="black", font=2) 
text(colMeans(bpctdown)-0.5 ,1, elementlabels, cex=1, pos=4, srt=90, font=1, xpd=NA)
mtext("Varied TRP element(s) \n(All others held fixed at intermediate level)", side=3, line=2)
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)

# calculate contribution of each element to total novel regimen impact
text(bpctdown[2,4], 140*final_pctdown[1,3,"0.5"], "Fraction of total variation in novel regimen's impact attributable to variation in specified characteristic", pos=1, xpd=NA)
text(bpctdown[2,2:7], 150*final_pctdown[1,3,"0.5"], paste0(round(((final_pctdown[,3,3]-final_pctdown[,1,3] )/(final_pctdown[1,3,3]-final_pctdown[1,1,3] ))[2:7]*100, digits=1), "%") , pos=1, xpd=NA)



# Same for DR regimen
## resize plot area to wide ##
novelwide <- allnovelwide[["DRDSTall"]]; drout <- alldrout[1:nrow(novelwide),]

outcome <- c("tbdeaths") 
y <- 10

r_pctdown <- array(0,dim=c( length(elementnames) , 3 , 5 )); 
dimnames(r_pctdown) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))

for (vary in elementnames) 
{ 
  r_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, y, vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  r_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, y,"allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  r_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, y, vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
}  


par(mar=c(3,5,8,3), mfrow=c(1,2))
ylow=250*r_pctdown[1,3,"0.5"]
bpctdown <- barplot(height = 100*aperm(r_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="Percent reduction in year 10 TB mortality with novel regimen,\n compared to projection under current standard of care", 
                    xlab="", las=2, cex.lab=1, 
                    ylim=c(1,-0.05)*ylow,
                    legend=c("minimal","intermediate","optimal"), args.legend=list(title="Level of varied element(s)",x=bpctdown[2,7], y=90*min(r_pctdown[1,3,"0.025"]),cex=0.8),
                    main="", 
                    col=cols, names.arg=rep("", length(elementnames)))
text(bpctdown+0.5, ylow/100, drlabels ,cex=0.9, pos=2, srt=90, col="black", font=2) 
text(colMeans(bpctdown)-0.5 ,-ylow/100, elementlabels, cex=0.8, pos=4, srt=90, font=1, xpd=NA)
mtext("Varied TRP element(s) \n(All others held fixed at intermediate level)", side=3, line=6)
arrows(bpctdown, aperm(100*r_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*r_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)

text(bpctdown[2,4], 7/8*ylow, "Fraction of total variation in novel regimen's impact\nattributable to variation in specified characteristic:", pos=1, xpd=NA)
text(bpctdown[2,2:7], 15/16*ylow, paste0(round(((r_pctdown[,3,3]-r_pctdown[,1,3] )/(r_pctdown[1,3,3]-r_pctdown[1,1,3] ))[2:7]*100, digits=1), "%") , pos=1, xpd=NA)


outcome <- c("rrdeaths") 

r_pctdown <- array(0,dim=c( length(elementnames) , 3 , 5 )); 
dimnames(r_pctdown) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))

for (vary in elementnames) 
{ 
  r_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  r_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  r_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
}  

ylow <- 130*r_pctdown[1,3,"0.5"]
bpctdown <- barplot(height = 100*aperm(r_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="Percent reduction in year 10 rifampin-resistant TB mortality", 
                    xlab="", las=2, cex.lab=1, 
                    ylim=c(1,-0.05)*ylow,
                    main="", 
                    col=cols, names.arg=rep("", length(elementnames)))
text(bpctdown+0.5, ylow/100, drlabels ,cex=0.9, pos=2, srt=90, col="black", font=2) 
text(colMeans(bpctdown)-0.5 ,-ylow/100, elementlabels, cex=0.8, pos=4, srt=90, font=1, xpd=NA)
mtext("Varied TRP element(s) \n(All others held fixed at intermediate level)", side=3, line=6)
arrows(bpctdown, aperm(100*r_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*r_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)

text(bpctdown[2,4], 7/8*ylow, "Fraction of total variation in novel regimen's impact\nattributable to variation in specified characteristic:", pos=1, xpd=NA)
text(bpctdown[2,2:7], 15/16*ylow, paste0(round(((r_pctdown[,3,3]-r_pctdown[,1,3] )/(r_pctdown[1,3,3]-r_pctdown[1,1,3] ))[2:7]*100, digits=1), "%") , pos=1, xpd=NA)

                    


############ duration

# sensitivity analysis: ltfu at 2 months
ltfu2mo <- read.csv("TRPwideoutput_DSDSTall_ltfu2mo.India_20160105.1.csv", header=TRUE)
drltfu <- screendrout("DRcalibration_ltfu2mo.India_20160105.1.csv", tolerance=1.5)
par(mfrow=c(1,1)); p<-plotpctdown(outcome="tbdeaths", scenario="DSDSTall", elements=elementnames, novelwide=ltfu2mo, drout =drltfu, barlabels=TRUE, elementlabs=TRUE, 
               main="With all losses to follow up occurring at 2 months\n(to maximize the impact of regimen duration)", mar=c(4,4,7,1))
text(p$positions[2,4], -34, "Fraction of total variation in novel regimen's impact\nattributable to variation in specified characteristic:", pos=1, xpd=NA)
text(p$positions[2,2:7], -37, paste0(round(((p$array[,3,3]-p$array[,1,3] )/(p$array[1,3,3]-p$array[1,1,3] ))[2:7]*100, digits=1), "%") , pos=1, xpd=NA)

# resource use: barplots as above but for outcomes of diagnoses, DSTs (rif and novel in same plot), and rxmonths (all 3 in same plot)
novelwide <- allnovelwide[["DSDSTall"]]; drout <- alldrout[1:nrow(novelwide),]

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


par(mar=c(4,6,3,2), mfcol=c(2,3), oma=c(2,1,1,1)) 
bres <- barplot(12*cbind(aperm(resource,c(1,3,2,4))["efficacy",1:3,,3], aperm(resource,c(1,3,2,4))["duration",1:3,,3], aperm(resource,c(1,3,2,4))["companion",1:3,,3]), beside = FALSE, 
                      space=c(0.75,0.25,0.25), cex.lab=1.2,
                    col=rev(blues), ylab="Patient-months of treatment\nin year 10, by regimen")
legend(x = -2*bres[1], y=max(12*rowSums(resource["efficacy",,1:3,3]))+120, xjust=0, yjust=0, fill=blues,
       c("novel regimen","second-line regimen","first-line regimen"), xpd=NA)
text(bres[c(2,5,8)], -5*12, elementlabels[2:4], cex=1, pos=1, xpd=NA)

bres <- barplot(12*cbind(aperm(cumresource,c(1,3,2))["efficacy",1:3,], aperm(cumresource,c(1,3,2))["duration",1:3,], aperm(cumresource,c(1,3,2))["companion",1:3,]), beside = FALSE, 
                space=c(0.75,0.25,0.25), cex.lab=1.2,
                col=rev(blues), ylab="Cumulative patient-months of treatment\nthrough year 10, by regimen")
text(bres[c(2,5,8)], -50*12, elementlabels[2:4], cex=1, pos=1, xpd=NA)


# diagnostic and other resources

tests <- array(0,dim=c( length(elementnames) , 3 , 3, 5 ));
dimnames(tests) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "test"=c("Novel regimen DSTs performed","Rifampin DSTs performed","Treatment courses initiated"),"q"=c(0.025,0.25,0.5,0.075,0.975))

for (vary in elementnames)  for (ntest in 1:3)
{ 
  test <- c("nDSTs","rDSTs","dxs")[ntest]
  tests[vary,1,ntest,] <- quantile(novelwide[ , paste0(test, "10", vary,"minimal")], c(0.025,0.25,0.5,0.075,0.975))
  tests[vary,2,ntest,] <- quantile(novelwide[ , paste0(test, "10allintermediate")], c(0.025,0.25,0.5,0.075,0.975))
  tests[vary,3,ntest,] <- quantile(novelwide[ , paste0(test, "10", vary,"optimal")], c(0.025,0.25,0.5,0.075,0.975))
}  


bup <- barplot(c(aperm(tests,c(1,3,2,4))["efficacy",3,,3], aperm(tests,c(1,3,2,4))["duration",3,,3], aperm(tests,c(1,3,2,4))["companion",3,,3]), beside = TRUE, 
                space=c(0.75,0.25,0.25), ylim=c(0,140), cex.lab=1.2,
                col=c("gray30","gray60","gray90"), ylab="Total treatment courses initiated, year 10", xlab="")
text(bup[c(2,5,8)], -12, elementlabels[2:4], cex=1, pos=1, xpd=NA)
# mtext("Varied TRP element (All others held fixed at intermediate level)", side=1, line=5, cex=0.7)
arrows(bup, aperm(tests, c(2,1,3,4))[,2:4,3,"0.025"], bup, aperm(tests, c(2,1,3,4))[,2:4,3,"0.975"], angle=90, code=3, length=0.05)


bup <- barplot(c(aperm(tests,c(1,3,2,4))["efficacy",1,,3], aperm(tests,c(1,3,2,4))["duration",1,,3], aperm(tests,c(1,3,2,4))["companion",1,,3]), beside = TRUE, 
               space=c(0.75,0.25,0.25), ylim=c(0,140), cex.lab=1.2,
               col=c("gray30","gray60","gray90"), ylab="Novel regimen DSTs performed, year 10", xlab="")
text(bup[c(2,5,8)], -12, elementlabels[2:4], cex=1, pos=1, xpd=NA)
# mtext("Varied TRP element (All others held fixed at intermediate level)", side=1, line=5, cex=0.7)
arrows(bup, aperm(tests, c(2,1,3,4))[,2:4,1,"0.025"], bup, aperm(tests, c(2,1,3,4))[,2:4,1,"0.975"], angle=90, code=3, length=0.05)

outcomes <- c("tbdeaths","inc","relapses") 

down <- array(0,dim=c( length(elementnames) , 3 , length(outcomes), 5 )); 
dimnames(down) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), outcome=outcomes, "q"=c(0.025,0.25,0.5,0.075,0.975))

for (vary in elementnames) for (nout in 1:length(outcomes))
{ 
  outcome <- outcomes[nout]
  down[vary,1,nout,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  down[vary,2,nout,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  down[vary,3,nout,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
}  

bdown <- barplot(100*c(down["efficacy",,1,3], down["duration",,1,3], down["companion",,1,3]), beside = TRUE, 
               space=c(0.75,0.25,0.25), ylim=c(-25,0), cex.lab=1.2,
               col=c("gray30","gray60","gray90"), ylab="% reduction in TB deaths, year 10", xlab="")
text(bup[c(2,5,8)], -28, elementlabels[2:4], cex=1, pos=1, xpd=NA)
# mtext("Varied TRP element (All others held fixed at intermediate level)", side=1, line=5, cex=0.7)
arrows(bdown, 100*c(t(down[2:4,,1,"0.025"])), bdown, 100*c(t(down[2:4,,1,"0.975"])), angle=90, code=3, length=0.05)

# bdown <- barplot(100*c(down["efficacy",,2,3], down["duration",,2,3], down["companion",,2,3]), beside = TRUE, 
#                  space=c(0.75,0.25,0.25), ylim=c(-25,0), cex.lab=1.2,
#                  col=c("gray30","gray60","gray90"), ylab="% reduction in TB incidence, year 10", xlab="")
# text(bup[c(2,5,8)], -28, elementlabels[2:4], cex=1, pos=1, xpd=NA)
# # mtext("Varied TRP element (All others held fixed at intermediate level)", side=1, line=5, cex=0.7)
# arrows(bdown, 100*c(t(down[2:4,,2,"0.025"])), bdown, 100*c(t(down[2:4,,2,"0.975"])), angle=90, code=3, length=0.05)
# 
tinc  <- array(0,dim=c( length(elementnames) , 3 , 5 )); 
dimnames(tinc) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))

for (vary in elementnames)
{ 
  tinc[vary,1,] <- quantile((novelwide[ , paste0("inc", "10", vary,"minimal")] + novelwide[ , paste0("relapses", "10", vary,"minimal")] - 
                               (drout[ , paste0("inc","10")]  + drout[ , paste0("relapses","10")] ))/
                              (drout[ , paste0("inc","10")] + drout[ , paste0("relapses","10")] ), c(0.025,0.25,0.5,0.075,0.975))
  tinc[vary,2,] <- quantile((novelwide[ , paste0("inc", "10allintermediate")] + novelwide[ , paste0("relapses", "10allintermediate")] - 
                               (drout[ , paste0("inc","10")] + drout[ , paste0("relapses","10")] ))/
                                   (drout[ , paste0("inc","10")] + drout[ , paste0("relapses","10")]), c(0.025,0.25,0.5,0.075,0.975))
  tinc[vary,3,] <- quantile((novelwide[ , paste0("inc", "10", vary,"optimal")] + novelwide[ , paste0("relapses", "10", vary,"optimal")] - 
                               (drout[ , paste0("inc","10")]  + drout[ , paste0("relapses","10")] ))/
                              (drout[ , paste0("inc","10")] + drout[ , paste0("relapses","10")] ), c(0.025,0.25,0.5,0.075,0.975))  
}

btinc <- barplot(100*c(tinc["efficacy",,3], tinc["duration",,3], tinc["companion",,3]), beside = TRUE, 
                 space=c(0.75,0.25,0.25), ylim=c(-15,0), cex.lab=1.2,
                 col=c("gray30","gray60","gray90"), ylab="% reduction in TB incidence, year 10", xlab="")
text(btinc[c(2,5,8)], -16.5, elementlabels[2:4], cex=1, pos=1, xpd=NA)
arrows(btinc, 100*c(t(tinc[2:4,,"0.025"])), btinc, 100*c(t(tinc[2:4,,"0.975"])), angle=90, code=3, length=0.05)


library("sensitivity")
prcc <- list()
prcc$duration <- pcc(X = novelwide[,40:73], y= (novelwide[,"tbdeaths10durationminimal"] - novelwide[,"tbdeaths10durationoptimal"] )/ (novelwide[,"tbdeaths10allminimal"] - novelwide[,"tbdeaths10alloptimal"] ), rank=TRUE)
cbind(rownames(prcc$duration$PRCC), prcc$duration$PRCC)[rev(order(abs(prcc$duration$PRCC))),]

# Novel DS regimen results are highly sensitive to the fraction of DR TB being correctly diagnosed, 
# and to outcomes of DR TB on both first and second line therapy, because we've assumed the novel regimen works for DS TB too even it's not the intended target. 

# Therefore
 ## Will redo India DS DSTall simulations assuming perfect rif DST (rDSTall), and maybe use that as primary dataset !!
par(mfrow=c(2,2))
novelwide <- read.csv("TRPwideoutput_DSDSTall_rDSTall.India_20160105.1.csv", header=TRUE)
drout <- alldrout[1:nrow(novelwide),]
plot(0:10, novelwide[1,paste0("rxtime_r",0:10,"allminimal")])


###############
#Companion drug resistance

novelwide <- allnovelwide[["DSDSTnone"]]; drout <- alldrout[1:nrow(novelwide),]

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

bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,c(4),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality", 
                    xlab="Baseline prevalence of\ncompanion-drug resistance", names.arg=c("10%","3%","0%"), cex.lab=1, main="With DST for\nnovel regimen",
                    ylim=c(-25,2), col=c("gray30","gray60","gray90"))
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,4,"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,4,"0.975"], angle=90, code=3, length=0.05)


bpctdown <- barplot(height = 100*aperm(final_pctdown2, c(2,1,3))[,c(4),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality", 
                    xlab="Baseline prevalence of\ncompanion-drug resistance", names.arg=c("10%","3%","0%"), cex.lab=1, main="Without DST for\nnovel regimen",
                    ylim=c(-25,2), col=c("gray30","gray60","gray90"))
arrows(bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,4,"0.025"], bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,4,"0.975"], angle=90, code=3, length=0.05)


# include DR

outcomes <- rep(c("tbdeaths", "panronsets", "rrdeaths", "panronsets"), each=2); 
scenarios <- c("DSDSTall", "DSDSTnone","DSDSTall", "DSDSTnone","DRDSTall","DRDSTnone","DRDSTall","DRDSTnone"); 
plotelements <- rep("companion",8)
par(mfcol=c(2,4))
for (version in 1:8) plotresult(outcomes[version], scenarios[version], plotelements[version])


## !! need to look at outcomes such as novelprev, rnovelprev once have new runs from MARCC (will replace current 1-2)
# change 2nd outcome above to novelprev or nprev?
# but note that there is the problem that c resistance isn't at true equilibrium, and without treatment to sustain companion resistance, cprev falls. But slowly. cnprev falls more rapidly 


#sensitivity
novelwide <- allnovelwide[["DSDSTall"]]; drout <- alldrout[1:nrow(novelwide),]
prcc$companion <- pcc(X = novelwide[,40:73], y= (novelwide[,"tbdeaths10companionminimal"] - novelwide[,"tbdeaths10companionoptimal"] )/ (novelwide[,"tbdeaths10allminimal"] - novelwide[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$companionres <- pcc(X = novelwide[,40:73], y= ( novelwide[,"novelprev10companionminimal"] / novelwide[,"novelprev10companionoptimal"] ), rank=TRUE)
cbind(rownames(prcc$companion$PRCC), prcc$companion$PRCC)[rev(order(abs(prcc$companion$PRCC))),]
cbind(rownames(prcc$companionres$PRCC), prcc$companionres$PRCC)[rev(order(abs(prcc$companionres$PRCC))),]
# summary: currently has the issue with dr's benefitting from "mis"treatment with the novel regimen, so will need to use rDSTall


novelwide <- allnovelwide[["DRDSTall"]]; drout <- alldrout[1:nrow(novelwide),]
prcc$companion_dr <- pcc(X = novelwide[,40:73], y= (novelwide[,"tbdeaths10companionminimal"] - novelwide[,"tbdeaths10companionoptimal"] )/ (novelwide[,"tbdeaths10allminimal"] - novelwide[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$companionres_dr <- pcc(X = novelwide[,40:73], y= ( novelwide[,"novelprev10companionminimal"] / novelwide[,"novelprev10companionoptimal"] ), rank=TRUE)
cbind(rownames(prcc$companion_dr$PRCC), prcc$companion_dr$PRCC)[rev(order(abs(prcc$companion_dr$PRCC))),]
cbind(rownames(prcc$companionres_dr$PRCC), prcc$companionres_dr$PRCC)[rev(order(abs(prcc$companionres_dr$PRCC))),]
# summary: correlated with parameters that increase the number of treatments (and to a lesser extent, that reduce the new dr transmissions to dilute the remaining c cases) per incident case




###########################
## Barrier to resistance ##

outcomes <- rep(c("tbdeaths", "panronsets", "rrdeaths", "panronsets"), each=2); 
scenarios <- c("DSDSTall", "DSDSTnone","DSDSTall", "DSDSTnone","DRDSTall","DRDSTnone","DRDSTall","DRDSTnone"); 
plotelements <- rep("barrier",8)
par(mfcol=c(2,4))
for (version in 1:8) plotresult(outcomes[version], scenarios[version], plotelements[version])


novelwide <- allnovelwide[["DSDSTall"]]; drout <- alldrout[1:nrow(novelwide),]
prcc$barrier <- pcc(X = novelwide[,40:73], y= (novelwide[,"tbdeaths10barrierminimal"] - novelwide[,"tbdeaths10barrieroptimal"] )/ (novelwide[,"tbdeaths10allminimal"] - novelwide[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$barrierres <- pcc(X = novelwide[,40:73], y= ( novelwide[,"novelprev10barrierminimal"] / novelwide[,"novelprev10barrieroptimal"] ), rank=TRUE)
cbind(rownames(prcc$barrier$PRCC), prcc$barrier$PRCC)[rev(order(abs(prcc$barrier$PRCC))),]
cbind(rownames(prcc$barrierres$PRCC), prcc$barrierres$PRCC)[rev(order(abs(prcc$barrierres$PRCC))),]


# synergy between barrier and companion? !!





############
# exclusions

# synergy with maximal efficacy
synergyheader <- c("inew", "ids","idr","targetprev","targetcoprev", "targetdr", "targetpt","DST", names(unlist(genericvalues)))
for (i in (2:length(elementnames))) synergyheader <- append(synergyheader,  #fixed this after the fact
                                                            paste0( rep(tallynames[c(1,8:23)], times=2*11), 
                                                                    rep( rep(0:10, each=length(tallynames)-6), 2),
                                                                    rep(elementnames[i], each=22*(length(tallynames)-6) ),
                                                                    rep( c("minimal","optimal"), each=11*(length(tallynames)-6) ) ) )

esynergy <- read.csv("efficacyhigh_DSDSTall_India_20160105.1.csv"); colnames(esynergy) <- synergyheader
drout <- alldrout[1:nrow(esynergy),]
outcome <- "tbdeaths"
edown <- array(0,dim=c( length(elementnames)-1 , 3, 5 ));
dimnames(edown) <- list("vary"=elementnames[-1], "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))

for (vary in elementnames[-1]) 
{ 
  edown[vary,1,] <- quantile((esynergy[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  edown[vary,2,] <- quantile((esynergy[ , paste0(outcome, "10efficacyoptimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  edown[vary,3,] <- quantile((esynergy[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
}  

par(mfrow=c(1,1), mar=c(9,3,4,1))
bpctdown <- barplot(height = 100*aperm(edown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    col=c("gray30","gray60","gray90"), ylim=100*min(edown[,,"0.5"])*c(1.5,-0.1),
                    names.arg=elementlabels[2:7], cex.names=0.7,
                    main="Exploring synergy with efficacy:\nNovel regimen impact with maximal % durably cured throughout\n(but for other characteristics still intermediate level)")
arrows(bpctdown, aperm(100*edown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*edown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)

text(mean(bpctdown), -32.5, "Fraction of total variation in novel regimen's impact attributable to variation in characteristic:
     \nWith novel regimen efficacy at maximal level (99% cure) throughout", pos=1, xpd=NA)
text(bpctdown[2,], -35.5, paste0(round((edown[,1,3]-edown[,3,3] )/(median(novelwide[,"tbdeaths10allminimal"]/drout[,"tbdeaths10"])-median(novelwide[,"tbdeaths10alloptimal"]/drout[,"tbdeaths10"]))*100, digits=1), "%") , pos=1, xpd=NA)
text(mean(bpctdown), -37, "With novel regimen efficacy at intermediate level (97% cure) throughout", pos=1, xpd=NA)
text(bpctdown[2,], -38, paste0(round(((final_pctdown[,3,3]-final_pctdown[,1,3] )/(final_pctdown[1,3,3]-final_pctdown[1,1,3] ))[2:7]*100, digits=1), "%") , pos=1, xpd=NA)


# synergy: impact of exclusion variation when efficacy intermediate vs when efficacy optimal: 
y <- c(final_pctdown["exclusions",1,3], final_pctdown["exclusions",3,3], edown["exclusions",1,3], edown["exclusions",3,3])*-100
b <- barplot(y, beside=TRUE, space=c(0.75,0,0.75,0), ylim=max(y)*c(1.5,-0.1))
text(x = b,y = y, round(y,1))
text(b, y/50, dslabels[c(16,18,16,18)]  ,cex=1, pos=2, srt=90, col="black", font=2)


bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    main=paste0("Novel ",substr(scenario,1,2)," TB regimen,\n",scenarionames[[scenario]]," DST"),
                    ylim=100*min(final_pctdown[,,"0.5"])* c(1.5,-0.2), ylab=outcomenames[[outcome]],
                    col=cols, xlab="")
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)
text(bpctdown+0.4, -0.1, dslabels[which(elements==elementnames) * 3 - (2:0)] ,cex=1, pos=2, srt=90, col="black", font=2)
mtext(paste0("Varying ",elementlabels[which(elements==elementnames)], "with fixed maximal efficacy"), side=1, cex=0.8, line=5)
  
