source("TPPmat.R")

currenttag <- "India_20160111"; tolerance <- 1.5; location=""
targetepi <- "India"

source("displayfunction.R")


###############################
par(mar=c(3,4,7,1))
bpctdown <- barplot(height = array(-c(1,3,5,2,3,4), dim=c(3,2)), beside = TRUE, 
                    ylab="Impact of novel regimen", ylim=c(-6,0), yaxt='n',
                    xlab="", las=1, cex.lab=1.5, 
                    legend=c("minimal","intermediate","optimal"), args.legend=list(title="Level of varied element(s)", x=8, y= 2.1, cex=1.2),
                    col=cols, names.arg=c("Over full range\n of possible novel regimens","Varying only\nspecified characteristic"), cex.names=1.1)
arrows(bpctdown[1,], -5.5, bpctdown[3,], -5.5, angle=90, code=3, length=0.1)


# Results overview #

novelwide <- allnovelwide[["DSDSTall_rDSTall."]]; drout <- droutds <- alldrDST[1:nrow(novelwide),]

outcome <- c("tbdeaths") #can set up loop over multiple outcomes

traj <- array(0, dim=c(11,3,5)); dimnames(traj) <- list("t"=0:10, "level"=levels, "q"=c(0.025,0.25,0.5,0.75,0.975))

for (t in 0:10) for (l in levels) traj[t+1,l,] <- quantile(novelwide[,colnames(novelwide)==paste0(outcome, t, "all",l)], c(0.025,0.25,0.5,0.75,0.975))

final_abs <- final_diff <- final_pct <- array(0,dim=c( length(elementnames) , 2 , 5 )); 
  final_down <- final_pctdown <- array(0,dim=c( length(elementnames) , 3 , 5 )); 
dimnames(final_abs) <- dimnames(final_diff) <- dimnames(final_pct) <-  list("vary"=elementnames, "level"=c("minimal", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
  dimnames(final_down) <- dimnames(final_pctdown) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))

for (vary in elementnames) 
{ 
  final_abs[vary,1,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"minimal")], c(0.025,0.25,0.5,0.75,0.975))
  final_abs[vary,2,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"optimal")], c(0.025,0.25,0.5,0.75,0.975))
  
  final_diff[vary,1,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"minimal")] - novelwide[ , paste0(outcome,"10allintermediate")], c(0.025,0.25,0.5,0.75,0.975))
  final_diff[vary,2,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"optimal")] - novelwide[ , paste0(outcome,"10allintermediate")], c(0.025,0.25,0.5,0.75,0.975))
  
  final_pct[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - novelwide[ , paste0(outcome,"10allintermediate")] )/
                                    novelwide[ , paste0(outcome,"10allintermediate")], c(0.025,0.25,0.5,0.75,0.975))
  final_pct[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - novelwide[ , paste0(outcome,"10allintermediate")] )/
                                    novelwide[ , paste0(outcome,"10allintermediate")], c(0.025,0.25,0.5,0.75,0.975))

  final_down[vary,1,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  final_down[vary,2,] <- quantile(novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  final_down[vary,3,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  
  final_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  final_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  final_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
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
                    ylab="Percent reduction (median [IQR]) in year 10 TB mortality with novel regimen,\n compared to projection under current standard of care", 
                    xlab="", las=2, cex.lab=1, 
                    ylim=c(150*final_pctdown[1,3,"0.5"], 8),
                    legend=c("minimal","intermediate","optimal"),
                    args.legend=list(title="Level of varied element(s)", x=25, y=-18,cex=0.8),
                    main="", space=c(0,0,0,1.5,0,0,rep(c(0.5,0,0),5)),
                    col=cols, names.arg=rep("", length(elementnames)))
text(bpctdown+0.4, -0.5, dslabels ,cex=0.9, pos=2, srt=90, col="black", font=2) 
text(colMeans(bpctdown)-0.5 ,1, elementlabels, cex=0.9, pos=4, srt=90, font=1, xpd=NA)
mtext("Varied TRP element(s)", side=3, line=2)
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.25"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.75"], angle=90, code=3, length=0.05)
# arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)
# calculate contribution of each element to total novel regimen impact
text(bpctdown[2,4], 140*final_pctdown[1,3,"0.5"], "Fraction of total variation in novel regimen's impact attributable to variation in specified characteristic", pos=1, xpd=NA)
text(bpctdown[3,2:7], 150*final_pctdown[1,3,"0.5"], paste0(round(((final_pctdown[,3,3]-final_pctdown[,1,3] )/(final_pctdown[1,3,3]-final_pctdown[1,1,3] ))[2:7]*100, digits=1), "%") , pos=1, xpd=NA)



# Same for DR regimen
## resize plot area to wide ##
novelwide3 <- novelwide <- allnovelwide[["DRDSTall_"]]; drout <- droutdr <- alldrout[1:nrow(novelwide),]

outcome <- c("tbdeaths") 

y <- 10
r_pctdown <- array(0,dim=c( length(elementnames) , 3 , 5 )); 
dimnames(r_pctdown) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))

for (vary in elementnames) 
{ 
  r_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, y, vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  r_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, y,"allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  r_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, y, vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}  


par(mar=c(3,5,8,3), mfrow=c(1,2))
ylow=300*r_pctdown[1,3,"0.5"]
bpctdown <- barplot(height = 100*aperm(r_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="% reduction, total TB mortality (median [IQR])", 
                    xlab="", las=2, cex.lab=1.3, 
                    ylim=c(1,-0.05)*ylow,
                    legend=c("minimal","intermediate","optimal"), args.legend=list(title="Level of varied element(s)",x=bpctdown[2,7], y=90*min(r_pctdown[1,3,"0.025"]),cex=0.8),
                    main="", space=c(0,0,0,1.5,0,0,rep(c(0.5,0,0),5)),
                    col=cols, names.arg=rep("", length(elementnames)))
text(bpctdown+0.5, ylow/100, drlabels ,cex=0.9, pos=2, srt=90, col="black", font=2) 
text(colMeans(bpctdown)-0.5 ,-ylow/100, elementlabels, cex=0.8, pos=4, srt=90, font=1, xpd=NA)
mtext("Varied TRP element(s)", side=3, line=6)
arrows(bpctdown, aperm(100*r_pctdown, c(2,1,3))[,,"0.25"], bpctdown, aperm(100*r_pctdown, c(2,1,3))[,,"0.75"], angle=90, code=3, length=0.05)


ylow <- 130*r_pctdown[1,3,"0.5"]
bpctdown <- barplot(height = 100*aperm(r_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="% reduction, DR (rifampin-resistant) TB mortality (median [IQR])", 
                    xlab="", las=2, cex.lab=1.2, 
                    ylim=c(1,-0.05)*ylow,
                    main="", space=c(0,0,0,1.5,0,0,rep(c(0.5,0,0),5)),
                    col=cols, names.arg=rep("", length(elementnames)))
text(bpctdown+0.5, ylow/100, drlabels ,cex=0.9, pos=2, srt=90, col="black", font=2) 
text(colMeans(bpctdown)-0.5 ,-ylow/100, elementlabels, cex=0.8, pos=4, srt=90, font=1, xpd=NA)
mtext("Varied TRP element(s)", side=3, line=6)
arrows(bpctdown, aperm(100*r_pctdown, c(2,1,3))[,,"0.25"], bpctdown, aperm(100*r_pctdown, c(2,1,3))[,,"0.75"], angle=90, code=3, length=0.05)

text(mean(bpctdown[2,4:5]), 7/8*ylow, "Fraction of total variation in novel regimen's impact\nattributable to variation in specified characteristic:", pos=1, xpd=NA)
text(bpctdown[3,2:7], 0.97*ylow, paste0(round((([,3,3]-r_pctdown[,1,3] )/(r_pctdown[1,3,3]-r_pctdown[1,1,3] ))[2:7]*100, digits=1), "%") , pos=1, xpd=NA)

                    

############ duration

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


par(mar=c(5,6,3,2), mfcol=c(1,1), oma=c(0,0,0,0)) 
bres <- barplot(12*cbind(aperm(resource,c(1,3,2,4))["efficacy",1:3,,3], aperm(resource,c(1,3,2,4))["duration",1:3,,3], aperm(resource,c(1,3,2,4))["companion",1:3,,3]), beside = FALSE, 
                      space=c(0.75,0.25,0.25), cex.lab=1.2, main="Treatment costs",
                    col=rev(blues), ylab="Patient-months of treatment\nin year 10, by regimen")
legend(x = bres[6], y=max(12*rowSums(resource["efficacy",,1:3,3]))+10, xjust=0, yjust=0, fill=blues,
       c("novel regimen","second-line regimen","first-line regimen"), xpd=NA)
text(bres[c(2,5,8)], -7*12, elementlabels[2:4], cex=0.9, pos=1, xpd=NA)

# bres <- barplot(12*cbind(aperm(cumresource,c(1,3,2))["efficacy",1:3,], aperm(cumresource,c(1,3,2))["duration",1:3,], aperm(cumresource,c(1,3,2))["companion",1:3,]), beside = FALSE, 
#                 space=c(0.75,0.25,0.25), cex.lab=1.2,
#                 col=rev(blues), ylab="Cumulative patient-months of treatment\nthrough year 10, by regimen")
# text(bres[c(2,5,8)], -50*12, elementlabels[2:4], cex=1, pos=1, xpd=NA)
# 

# diagnostic and other resources

tests <- array(0,dim=c( length(elementnames) , 3 , 3, 5 ));
dimnames(tests) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "test"=c("Novel regimen DSTs performed","Rifampin DSTs performed","Treatment courses initiated"),"q"=c(0.025,0.25,0.5,0.75,0.975))

for (vary in elementnames)  for (ntest in 1:3)
{ 
  test <- c("nDSTs","rDSTs","dxs")[ntest]
  tests[vary,1,ntest,] <- quantile((1-novelwide$initialloss_s)*novelwide[ , paste0(test, "10", vary,"minimal")], c(0.025,0.25,0.5,0.75,0.975))
  tests[vary,2,ntest,] <- quantile((1-novelwide$initialloss_s)*novelwide[ , paste0(test, "10allintermediate")], c(0.025,0.25,0.5,0.75,0.975))
  tests[vary,3,ntest,] <- quantile((1-novelwide$initialloss_s)*novelwide[ , paste0(test, "10", vary,"optimal")], c(0.025,0.25,0.5,0.75,0.975))
}  

par(mar=c(5,5.5,1,0),oma=c(0.5,0.5,0.5,0.5), mfrow=c(1,1))
bup <- barplot(c(aperm(tests,c(1,3,2,4))["efficacy",3,,3], aperm(tests,c(1,3,2,4))["duration",3,,3], aperm(tests,c(1,3,2,4))["companion",3,,3]), beside = TRUE, 
                space=c(0.75,0.25,0.25), ylim=c(0,100), cex.lab=1.2,
                col=c("gray30","gray60","gray90"), ylab="Total treatment courses\ninitiated, year 10", xlab="")
text(bup[c(2,5,8)], -20, elementlabels[2:4], cex=1, pos=1, xpd=NA)
arrows(bup, aperm(tests, c(2,1,3,4))[,2:4,3,"0.25"], bup, aperm(tests, c(2,1,3,4))[,2:4,3,"0.75"], angle=90, code=3, length=0.05)

bup <- barplot(c(aperm(tests,c(1,3,2,4))["efficacy",1,,3], aperm(tests,c(1,3,2,4))["duration",1,,3], aperm(tests,c(1,3,2,4))["companion",1,,3]), beside = TRUE, 
               space=c(0.75,0.25,0.25), ylim=c(0,70), cex.lab=1.2,
               col=cols, ylab="Novel regimen DSTs\nperformed, year 10", xlab="")
text(bup[c(2,5,8)], -15, elementlabels[2:4], cex=0.9, pos=1, xpd=NA)
arrows(bup, aperm(tests, c(2,1,3,4))[,2:4,1,"0.25"], bup, aperm(tests, c(2,1,3,4))[,2:4,1,"0.75"], angle=90, code=3, length=0.05)

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

# bdown <- barplot(100*c(down["efficacy",,2,3], down["duration",,2,3], down["companion",,2,3]), beside = TRUE, 
#                  space=c(0.75,0.25,0.25), ylim=c(-25,0), cex.lab=1.2,
#                  col=c("gray30","gray60","gray90"), ylab="% reduction in TB incidence, year 10", xlab="")
# text(bup[c(2,5,8)], -28, elementlabels[2:4], cex=1, pos=1, xpd=NA)
# # mtext("Varied TRP element (All others held fixed at intermediate level)", side=1, line=5, cex=0.7)
# arrows(bdown, 100*c(t(down[2:4,,2,"0.025"])), bdown, 100*c(t(down[2:4,,2,"0.975"])), angle=90, code=3, length=0.05)
# 
tinc  <- array(0,dim=c( length(elementnames) , 3 , 5 )); 
dimnames(tinc) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))

for (vary in elementnames)
{ 
  tinc[vary,1,] <- quantile((novelwide[ , paste0("inc", "10", vary,"minimal")] + novelwide[ , paste0("relapses", "10", vary,"minimal")] - 
                               (drout[ , paste0("inc","10")]  + drout[ , paste0("relapses","10")] ))/
                              (drout[ , paste0("inc","10")] + drout[ , paste0("relapses","10")] ), c(0.025,0.25,0.5,0.75,0.975))
  tinc[vary,2,] <- quantile((novelwide[ , paste0("inc", "10allintermediate")] + novelwide[ , paste0("relapses", "10allintermediate")] - 
                               (drout[ , paste0("inc","10")] + drout[ , paste0("relapses","10")] ))/
                                   (drout[ , paste0("inc","10")] + drout[ , paste0("relapses","10")]), c(0.025,0.25,0.5,0.75,0.975))
  tinc[vary,3,] <- quantile((novelwide[ , paste0("inc", "10", vary,"optimal")] + novelwide[ , paste0("relapses", "10", vary,"optimal")] - 
                               (drout[ , paste0("inc","10")]  + drout[ , paste0("relapses","10")] ))/
                              (drout[ , paste0("inc","10")] + drout[ , paste0("relapses","10")] ), c(0.025,0.25,0.5,0.75,0.975))  
}

btinc <- barplot(100*c(tinc["efficacy",,3], tinc["duration",,3], tinc["companion",,3]), beside = TRUE, 
                 space=c(0.75,0.25,0.25), ylim=c(-15,2), cex.lab=1.2,
                 col=cols, ylab="% reduction in TB incidence, year 10", xlab="")
text(btinc[c(2,5,8)], -16.5, elementlabels[2:4], cex=1, pos=1, xpd=NA)
arrows(btinc, 100*c(t(tinc[2:4,,"0.025"])), btinc, 100*c(t(tinc[2:4,,"0.975"])), angle=90, code=3, length=0.05)



# sensitivity analysis: ltfu at 2 months
ltfu2mo <- read.csv("TRPwideoutput_DSDSTall_ltfu2mo.rDSTall.India_20160111.1.csv", header=TRUE)
drltfu <- screendrout("DRcalibration_ltfu2mo.rDSTall.India_20160111.1.csv", tolerance=1.5)[1:nrow(ltfu2mo),]
par(mfrow=c(1,1), oma=c(0,0,2.5,0)); p<-plotpctdown(outcome="tbdeaths", scenario="DSDSTall_rDSTall.", elements=elementnames, novelwide=ltfu2mo, drout =drltfu, barlabels=TRUE, elementlabs=TRUE, 
                                  mar=c(0,4,5,1), main="")
mtext("With all losses to follow up occurring at 2 months\n(to maximize the impact of regimen duration)",side=3, outer=TRUE, cex=1.2, font=2, line=0) 


efracs <- (p$array[2:7,3,3]-p$array[2:7,1,3] )/sum(p$array[2:7,3,3]-p$array[2:7,1,3] )
xp <- numeric(6); for (i in 1:6) xp[i] <- sum(efracs[1:i])
xp2 <- c(0,xp[-6]) + (xp - c(0,xp[-6]))/2

ffracs <- (final_pctdown[2:7,3,3]-final_pctdown[2:7,1,3] )/sum(final_pctdown[2:7,3,3]-final_pctdown[2:7,1,3] )
fp <- numeric(6); for (i in 1:6) fp[i] <- sum(ffracs[1:i])
fp2 <- c(0,fp[-6]) + (fp - c(0,fp[-6]))/2

par(mfrow=c(1,1), mar=c(3,9,3,1), oma=c(0,0,0,0)) 
b <- barplot(array(c(efracs, ffracs), dim=c(6,2)), horiz = TRUE, beside=FALSE, las=2, 
              names.arg=c("With evenly-\ndistributed loss\nto follow up","With all losses to\nfollow up occurring\nat 2 months"), font=2,
             legend.text=elementnames[-1], col=rainbow(6), args.legend=list(x="left"))
mtext("Fraction of total variation in novel regimen's impact attributable to\nvariation in specified characteristic", side=3, font=2, cex=1.2)
text(rbind(xp2,fp2)-0.005,b+0.1*(rep(c(0,0,1,1),3)),round(rbind(efracs,ffracs),2))


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
b <- barplot(prcc$durationres$PRCC$original[display], horiz=TRUE, beside=TRUE, space=c(0.5,0),
             names.arg=longparamnames[display], 
             axes=FALSE,axisnames=FALSE, xlim=c(-1,1),
             las=1, cex.names=0.8, cex.axis=0.8,cex.lab=0.8, main="Novel DR regimen", cex.main=1 )
axis(1, at=0.2*c(-4:4),cex.axis=0.8,cex.lab=0.9)
mtext("Partial rank correlation of model parameters with the impact of regimen duration", side=3,cex=1.2,outer=TRUE)

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


novelwide <- novelwide2 <- allnovelwide[["DSDSTnone_rDSTall."]]; drout <- alldrDST[1:nrow(novelwide),]

final_pctdown2 <- array(0,dim=c( length(elementnames) , 3 , 5 ));
dimnames(final_pctdown2) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))

# and same for DR regimen:# 
novelwide <- novelwide2 <- allnovelwide[["DRDSTnone_"]]; drout <- alldrout[1:nrow(novelwide),]

outcome <- "rrdeaths"

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

bpctdown <- barplot(height = 100*aperm(r_pctdown, c(2,1,3))[,c(4),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 rifampin-resistant\nTB mortality (median [IQR])", 
                    xlab="Baseline prevalence of\ncompanion-drug resistance", names.arg=c("15%","5%","0%"), cex.lab=1, 
                    main="With DST for\nnovel DR regimen",
                    ylim=c(-40,2), col=cols)
arrows(bpctdown, aperm(100*r_pctdown, c(2,1,3))[,4,"0.25"], bpctdown, aperm(100*r_pctdown, c(2,1,3))[,4,"0.75"], angle=90, code=3, length=0.05)

bpctdown <- barplot(height = 100*aperm(final_pctdown2, c(2,1,3))[,c(4),"0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 rifampin-resistant\nTB mortality(median [IQR])", 
                    xlab="Baseline prevalence of\ncompanion-drug resistance", names.arg=c("15%","5%","0%"), cex.lab=1, 
                    main="Without DST for\nnovel DR regimen",
                    ylim=c(-40,2), col=cols)
arrows(bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,4,"0.25"], bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,4,"0.75"], angle=90, code=3, length=0.05)


# include DR

outcomes <- rep(c("tbdeaths", "novelprev", "rrdeaths", "novelprev"), each=2); 
scenarios <- c("DSDSTall_rDSTall.", "DSDSTnone_rDSTall.","DSDSTall_rDSTall.", "DSDSTnone_rDSTall.","DRDSTall_","DRDSTnone_","DRDSTall_","DRDSTnone_"); 
plotelements <- rep("companion",8)
par(mfcol=c(2,4))
for (version in 1:8) plotresult(outcomes[version], scenarios[version], plotelements[version])


## !! need to look at outcomes such as novelprev, rnovelprev 
# and trends over time in each??
# but note that there is the problem that c resistance isn't at true equilibrium, and without treatment to sustain companion resistance, cprev falls. But slowly. cnprev falls more rapidly 


#sensitivity
novelwide <- allnovelwide[["DSDSTall_rDSTall."]]; drout <- alldrDST[1:nrow(novelwide),]
prcc$companion <- pcc(X = novelwide[,40:73], y= (novelwide[,"tbdeaths10companionminimal"] - novelwide[,"tbdeaths10companionoptimal"] )/ (novelwide[,"tbdeaths10allminimal"] - novelwide[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$companionres <- pcc(X = novelwide[,40:73], y= ( novelwide[,"cprev10companionminimal"] / novelwide[,"cprev10companionoptimal"] ), rank=TRUE)
cbind(rownames(prcc$companion$PRCC), prcc$companion$PRCC)[rev(order(abs(prcc$companion$PRCC))),]
cbind(rownames(prcc$companionres$PRCC), prcc$companionres$PRCC)[rev(order(abs(prcc$companionres$PRCC))),]
# summary: 

novelwide <- allnovelwide[["DRDSTall_"]]; drout <- alldrout[1:nrow(novelwide),]
prcc$companion_dr <- pcc(X = novelwide[,40:73], y= (novelwide[,"tbdeaths10companionminimal"] - novelwide[,"tbdeaths10companionoptimal"] )/ (novelwide[,"tbdeaths10allminimal"] - novelwide[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$companionres_dr <- pcc(X = novelwide[,40:73], y= ( novelwide[,"novelprev10companionminimal"] / novelwide[,"novelprev10companionoptimal"] ), rank=TRUE)
cbind(rownames(prcc$companion_dr$PRCC), prcc$companion_dr$PRCC)[rev(order(abs(prcc$companion_dr$PRCC))),]
cbind(rownames(prcc$companionres_dr$PRCC), prcc$companionres_dr$PRCC)[rev(order(abs(prcc$companionres_dr$PRCC))),]
# summary: correlated with parameters that increase the number of treatments (and to a lesser extent, that reduce the new dr transmissions to dilute the remaining c cases) per incident case

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
b <- barplot(prcc$durationres$PRCC$original[display], horiz=TRUE, beside=TRUE, space=c(0.5,0),
             names.arg=longparamnames[display], 
             axes=FALSE,axisnames=FALSE, xlim=c(-1,1),
             las=1, cex.names=0.8, cex.axis=0.8,cex.lab=0.8, main="Novel DR regimen", cex.main=1 )
axis(1, at=0.2*c(-4:4),cex.axis=0.8,cex.lab=0.9)
mtext("Partial rank correlation of model parameters with the impact of regimen duration", side=3,cex=1.2,outer=TRUE)



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
dimnames(edown) <- list("vary"=elementnames[-1], "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))

for (vary in elementnames[-1]) 
{ 
  edown[vary,1,] <- quantile((esynergy[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  edown[vary,2,] <- quantile((esynergy[ , paste0(outcome, "10efficacyoptimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  edown[vary,3,] <- quantile((esynergy[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
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
  
