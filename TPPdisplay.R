source("TPPmat.R")

currenttag <- "India_20160313p"; tolerance <- 1.5; location=""
targetepi <- "India"

source("loaddata.R")


###############################
par(mar=c(3,6,7,1), mfrow=c(1,1),oma=c(0,0,0,0))
bpctdown <- barplot(height = array(-c(1,3,5,2,3,4), dim=c(3,2)), beside = TRUE, 
                    ylab="Impact of novel regimen\n(mortality reduction)", ylim=c(-6,0), yaxt='n',
                    xlab="", las=1, cex.lab=1.5, 
                    legend=c("minimal","intermediate","optimal"), args.legend=list(title="Level of varied element(s)", x=8, y= 2.2, cex=1.2),
                    col=cols, names.arg=c("Over full range\n of possible novel regimens","Varying only\nspecified characteristic"), cex.names=1.1)
arrows(bpctdown[1,], -5.5, bpctdown[3,], -5.5, angle=90, code=3, length=0.1)



# Results overview #

novelwide <- allnovelwide[["DSDSTall_rDSTall."]]; drout <- droutds <- alldrDST #[1:nrow(novelwide),]
outcome <- "tbdeaths"

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


final_pctdown <- array(0,dim=c( length(dselementnames) , 3 , 5 )); 
dimnames(final_pctdown) <- list("vary"=dselementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in dselementnames) 
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
                    space=c(0,0,0,1.5,0,0,rep(c(0.5,0,0),6)),
                    col=cols, names.arg=rep("", length(dselementnames)))
axis(2, at=seq(-25,0,by=5), labels=paste0(seq(-25,0,by=5),"%"),las=2,cex.axis=0.8)
mtext("Reduction in year 10 TB mortality with a novel DS-TB regimen,\n compared to projection under current standard of care", cex=1.2, line=5, side=3)
mtext("Varied TRP element(s)", side=3, line=3)
# text((bpctdown+0.4)[,2], -0.5, dslabels[4:6] ,cex=0.9, pos=2, srt=90, col="black", font=2) 
text(colMeans(bpctdown) ,0.5, dselementlabels, cex=0.8, pos=3, srt=0, font=1, xpd=NA)
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.25"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.75"], angle=90, code=3, length=0.05, xpd=NA)
# arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)
# calculate contribution of each element to total novel regimen impact

# # box version
# int <- final_pctdown["all",2,"0.5"]; min <- final_pctdown[,1,"0.5"] - int; opt <- final_pctdown[,3,"0.5"] - int
# par(mar=c(3,3,3,12))
# boxes <- barplot(100*(min-int), offset=100*int, ylim=100*1.2*c(min(opt),max(min)),
#                 ylab="% Change in TB mortality",
#                 col=cols[1], names=dselementnames)
# barplot(100*(opt-int), offset=100*int, col=cols[3], add=TRUE)
# abline(h=0,lty=2); text(boxes[8]+1,0.5,"Baseline projection under\ncontinuing current care", xpd=NA, pos=4)
# abline(h=100*int,lty=2); text(boxes[8]+1,100*int+0.5,"Estimated impact of novel\nregimen meeting all\nintermediate targets", xpd=NA, pos=4)
# legend("bottomright", c("Characteristic(s) kept at minimal level", "Characteristic(s) improved to optimal level"),
#        fill=cols[c(1,3)], xjust=0)

## incidence outcome

outcome <- c("inc") #can set up loop over multiple outcomes

incdown <- array(0,dim=c( length(dselementnames) , 3 , 5 )); 
dimnames(incdown) <- list("vary"=dselementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in dselementnames) 
{ incdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  incdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  incdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                       drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}  

par(mar=c(2,4,8,1), mfrow=c(1,1), oma=c(0,0,0,0))
bpctdown <- barplot(height = 100*aperm(incdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="% reduction (median [IQR])",
                    xlab="", las=2, cex.lab=1, 
                    ylim=c(150*incdown[1,3,"0.5"], 2), yaxt='n',
                    legend=c("minimal","intermediate","optimal"),
                    args.legend=list(title="Level of varied element(s)", x="bottom",cex=0.9),
                    space=c(0,0,0,1.5,0,0,rep(c(0.5,0,0),6)),
                    col=cols, names.arg=rep("", length(dselementnames)))
axis(2, at=seq(-12,0,by=2), labels=paste0(seq(-12,0,by=2),"%"),las=2,cex.axis=0.8)
mtext(expression(paste("Reduction in year 10 TB ",bold("incidence")," with a novel DS-TB regimen,")), cex=1.2, line=5, side=3)
mtext("compared to projection under current standard of care", cex=1.2, line=4, side=3)
mtext("Varied TRP element(s)", side=3, line=2)
text(colMeans(bpctdown) ,0.5, dselementlabels, cex=0.8, pos=3, srt=0, font=1, xpd=NA)
arrows(bpctdown, aperm(100*incdown, c(2,1,3))[,,"0.25"], bpctdown, aperm(100*incdown, c(2,1,3))[,,"0.75"], angle=90, code=3, length=0.05, xpd=NA)


# DS combos

novelwide <- allnovelwide$DSDSTall_rDSTall.
dscombos <- numeric(0)
i <- 1; while(file.exists(paste0(location,"TRPcombos_DSDSTall_rDSTall.",currenttag,".",i,".csv")))
{dscombos <- rbind(dscombos, read.csv(paste0(location,"TRPcombos_DSDSTall_rDSTall.",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0read.csv(paste0())
outcome <- "tbdeaths"
drDST <- alldrDST

scenarionames <- c("em_dm_tm_xm_sm_","em_dm_tm_xm_si_","em_dm_tm_xm_so_", 
                   "em_do_tm_xm_sm_", "em_do_tm_xm_si_", "em_do_tm_xm_so_", 
                   "em_dm_to_xm_sm_", "em_dm_to_xm_si_", "em_dm_to_xm_so_", 
                   "em_dm_tm_xo_sm_", "em_dm_tm_xo_si_", "em_dm_tm_xo_so_",
                   
                   "ei_dm_tm_xm_sm_", "ei_dm_tm_xm_si_", "ei_dm_tm_xm_so_", 
                   "ei_do_tm_xm_sm_","ei_do_tm_xm_si_", "ei_do_tm_xm_so_", 
                   "ei_dm_to_xm_sm_", "ei_dm_to_xm_si_", "ei_dm_to_xm_so_", 
                   "ei_dm_tm_xo_sm_","ei_dm_tm_xo_si_", "ei_dm_tm_xo_so_",
                   
                   "eo_dm_tm_xm_sm_","eo_dm_tm_xm_si_", "eo_dm_tm_xm_so_",
                   "eo_do_tm_xm_sm_","eo_do_tm_xm_si_", "eo_do_tm_xm_so_", 
                   "eo_dm_to_xm_sm_","eo_dm_to_xm_si_", "eo_dm_to_xm_so_", 
                   "eo_dm_tm_xo_sm_","eo_dm_tm_xo_si_", "eo_dm_tm_xo_so_")

a <- numeric(0); for (i in 1:length(grep("^tbdeaths",colnames(dscombos)))) a[i] <- median((dscombos[,grep("^tbdeaths",colnames(dscombos))[i]]- drDST$tbdeaths10 )/drDST$tbdeaths10)
a[length(a)+1] <- median(((allnovelwide$DSDSTall_rDSTall.$tbdeaths10alloptimal - alldrDST$tbdeaths10 )/alldrDST$tbdeaths10))
names(a) <- c(colnames(dscombos)[grep("^tbdeaths",colnames(dscombos))], "all optimal")
par(mar=c(3,5,3,3)); 
bcols <- c("gray",topo.colors(3))
b<-barplot(a[c(c(1,10,7,4,12,9,6)+12*rep(c(0,1,2), each=7))], las=2, ylab="TB mortality reduction", space = c(rep(0.1,7),1,rep(0.1,6),1,rep(0.1,6)),
           col=c( rep(c( bcols[1], rep(bcols[2:4], 2)), times=3)), angle=20, density=c( rep( c(rep(30,4), rep(100,3)), times=3), 100) ,
           ylim=c(-0.16,0.02), xaxt='n', cex.lab=1.2,
           main="Impact of select characteristic combinations,\nwith and without associated increased novel regimen uptake")
text(x=b[c(4,11,18)], y=0.01, c("Current efficacy", "Intermediate efficacy", "Optimal efficacy"), cex=1.2)
legend("bottomleft", c("Baseline ~ standard of care", "Increased eligibility", "Improved adherence", "Shortened duration"),
       fill=bcols[1:4], xpd=NA, bty="n", cex=1.2)
legend("bottomright", c("With limited (50%) uptake of novel DS regimen", "With universal uptake of novel DS regimen"), 
       angle=20, density = c(30,100), fill = c("black","black"), bty='n', cex=1.2)




# Minimum and maximal baselines
## novelwide's allminimal and alloptimal include variation in scaleup, so "allbut" ends up making extra difference because it also includes loss of optimal scaleup. so use allminopt runs instead.

novelwide <- allnovelwide[["DSDSTall_rDSTall."]]; drout <- droutds
outcome <- "tbdeaths"

###!! or should I benchmark to drout[,paste0(outcome,"10")] instead of allminopt_ds[ , paste0(outcome, "allmin")]? I think here benchmarking to an all-minimal regimen makes sense.

dsmin <- array(0,dim=c(6,5)); dimnames(dsmin) <- list("element"=dselementnames[-c(1,8)], "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in dselementnames[2:7])
{
  dsmin[vary,] <- quantile( 
     (allminopt_ds[ , paste0(outcome, "10allmin")] - only_ds[ , paste0(outcome, "10only",vary)] ) /
      (allminopt_ds[ , paste0(outcome, "10allmin")] - allminopt_ds[ , paste0(outcome, "10allopt")]) , c(0.025,0.25,0.5,0.75,0.975) )
}

par(mar=c(1,1,3,1), mfrow=c(1,1)) 
b <- barplot(dsmin[,3], horiz = TRUE, beside=TRUE, las=2, font=2, xaxt='n',
             legend.text=shortelementlabels[-c(1,8,9)], col=rainbow(6), args.legend=list(x="topright", cex=1), xlim=c(-0.2,1))
mtext("Fraction of total mortality benefit of all-optimal novel DS regimen (median [95% UR]) that is
      attainable by improving only a single characteristic from minimal to optimal level", side=3, line=1, font=2, cex=1, xpd=NA)
text(dsmin[,5]+0.05, b, paste0(round(dsmin[,3]*100, 1),"%"))
arrows(dsmin[,"0.025"], b, dsmin[,"0.975"], b, angle=90, code=3, length=0.05)
text(0.45, b[2], "'All minimal' baseline here is a regimen with 
      ~ standard-of-care efficacy, duration, and tolerability,
      that excludes 10% due to resistance, that excludes
      10% of patients for medical reasons,
      that reaches only 50% of eligible patients,
      and that has low barrier to resistance.", pos=4)

novelwide <- allnovelwide[["DSDSTall_rDSTall."]]; drout <- droutds
outcome <- "tbdeaths"
intdsmin <- array(0,dim=c(6,5)); dimnames(intdsmin) <- list("element"=dselementnames[-c(1,8)], "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in dselementnames[2:7])
{
  intdsmin[vary,] <- quantile( 
    (allminopt_ds[ , paste0(outcome, "10allmin")] - intonly_ds[ , paste0(outcome, "10only",vary)] ) /
      (allminopt_ds[ , paste0(outcome, "10allmin")] - allminopt_ds[ , paste0(outcome, "10allopt")]) , c(0.025,0.25,0.5,0.75,0.975) )
}  

novelwide <- allnovelwide$DSDSTall_rDSTall.
drout <- alldrDST
outcome <- "tbdeaths"

dsmax <- array(0,dim=c(6,5)); dimnames(dsmax) <- list("element"=dselementnames[-c(1,8)], "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in dselementnames[2:7])
{
  dsmax[vary,] <- quantile( 
      1- ( allbut_ds[ , paste0(outcome, "10allbut",vary)] -  drout[ , paste0(outcome,"10")]) /
          (allminopt_ds[ , paste0(outcome, "10allopt")] - drout[ , paste0(outcome,"10")]) , c(0.025,0.25,0.5,0.75,0.975) )
}

par(mar=c(1,1,3,1), mfrow=c(1,1)) 
b <- barplot(dsmax[,3], horiz = TRUE, beside=TRUE, las=2, font=2, xaxt='n', xlim=c(0,1),
             legend.text=shortelementlabels[-c(1,8,9)], col=rainbow(6), args.legend=list(x="topright", cex=1))#1.5, y=8, cex=0.8))
mtext("Fraction of an all-optimal novel DS regimen's mortality impact (median [95% UR])\nthat is lost when only one characteristic is reduced to its minimal value", side=3, line=1, font=2, cex=1.1, xpd=NA)
arrows(dsmax[,"0.025"], b, dsmax[,"0.975"], b, angle=90, code=3, length=0.05)
text(dsmax[,5]+0.05, b, paste0(round(dsmax[,3]*100, 1),"%"))

novelwide <- allnovelwide$DSDSTall_rDSTall.
drout <- droutds
outcome <- "tbdeaths"

intdsmax <- array(0,dim=c(6,5)); dimnames(intdsmax) <- list("element"=dselementnames[-c(1,8)], "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in dselementnames[2:7])
{
  intdsmax[vary,] <- quantile( 
    1- ( intallbut_ds[ , paste0(outcome, "10allbut",vary)] -  drout[ , paste0(outcome,"10")]) /
      (allminopt_ds[ , paste0(outcome, "10allopt")] - drout[ , paste0(outcome,"10")]) , c(0.025,0.25,0.5,0.75,0.975) )
}

par(mar=c(1,1,3,1), mfrow=c(1,1)) 
b <- barplot(intdsmax[,3], horiz = TRUE, beside=TRUE, las=2, font=2, xaxt='n', xlim=c(0,1),
             legend.text=shortelementlabels[-c(1,8,9)], col=rainbow(6), args.legend=list(x="topright", cex=1))#1.5, y=8, cex=0.8))
mtext("Fraction of an all-optimal novel DS regimen's mortality impact (median [95% UR])\nthat is lost when only one characteristic is reduced to its intermediate value", side=3, line=1, font=2, cex=1.1, xpd=NA)
arrows(intdsmax[,"0.025"], b, intdsmax[,"0.975"], b, angle=90, code=3, length=0.05)
text(intdsmax[,5]+0.05, b, paste0(round(intdsmax[,3]*100, 1),"%"))

alldsmin <- array(0,dim=c(12,6)); for (i in 1:6) alldsmin[((1:2)+(2*(i-1))),i] <- rbind(intdsmin[,3], dsmin[,3]-intdsmin[,3])[,i]
alldsmax <- array(0,dim=c(12,6)); for (i in 1:6) alldsmax[((1:2)+(2*(i-1))),i] <- rbind(intdsmax[,3], dsmax[,3]-intdsmax[,3])[,i]

par(mar=c(1,1,3,1), mfrow=c(1,1)) 
b <- barplot(alldsmax, horiz = TRUE, beside=FALSE, las=2, font=2, xaxt='n', xlim=c(0,1),
             col=rep(rainbow(6),each=2), density=c(200,25))
legend(x="topright", cex=1, legend=c(shortelementlabels[-c(1,8,9)],"Reduced to intermediate","Reduced to minimal"), fill=c(rainbow(6),"black","black"), density=c(rep(200,7),25))
mtext("Fraction of an all-optimal novel DS regimen's mortality impact (median [95% UR])\nthat is lost when only one characteristic is reduced to the intermediate or minimal target", side=3, line=1, font=2, cex=1.1, xpd=NA)
# arrows(intdsmax[,"0.025"], b, intdsmax[,"0.975"], b, angle=90, code=3, length=0.05)
arrows(dsmax[,"0.025"], b, dsmax[,"0.975"], b, angle=90, code=3, length=0.05)
text(dsmax[,5]+0.05, b, paste0(round(dsmax[,3]*100, 1),"%"))

# Same for DR regimen

## resize plot area to wide ##
novelwide <- novelwide3 <- allnovelwide[["DRDSTall_"]]; drout <- droutdr <- alldrout#[1:nrow(novelwide),]

outcome <- c("tbdeaths") 
r_pctdown <- array(0,dim=c( length(drelementnames) , 3 , 5 )); y <- 10
dimnames(r_pctdown) <- list("vary"=drelementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in drelementnames) 
{ r_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, y, vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  r_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, y,"allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  r_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, y, vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                   drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}  

outcome <- c("rrdeaths"); y <- 10
rr_pctdown <- array(0,dim=c( length(drelementnames) , 3 , 5 )); 
dimnames(rr_pctdown) <- list("vary"=drelementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in drelementnames) 
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
                    space=c(0,0,0,1.5,0,0,rep(c(0.5,0,0),6)),
                    col=cols, names.arg=rep("", length(drelementnames)))
axis(2, at=c(-5,0), labels=paste0(c(5,0),"%"),las=2,cex.axis=0.8)
text(colMeans(bpctdown) ,ylow, shortelementlabels[c(1:7,9)], cex=0.8, pos=4, srt=90, font=2, xpd=NA)
mtext("Varied TRP element(s)", side=1, line=1, xpd=NA, cex=0.8)
arrows(bpctdown, aperm(100*r_pctdown, c(2,1,3))[,,"0.25"], bpctdown, aperm(100*r_pctdown, c(2,1,3))[,,"0.75"], angle=90, code=3, length=0.05, xpd=NA)

ylow <- 130*rr_pctdown[1,3,"0.5"]
bpctdown <- barplot(height = 100*aperm(rr_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="% reduction, DR-TB mortality (median [IQR])", 
                    main="", xlab="", las=2, cex.lab=1, 
                    ylim=c(1,0)*ylow, yaxt='n',
                    space=c(0,0,0,1.5,0,0,rep(c(0.5,0,0),6)),
                    col=cols, names.arg=rep("", length(drelementnames)))
axis(2, at=-seq(0,-ylow,by=5), labels=paste0(seq(0,-ylow,by=5),"%"),las=2,cex.axis=0.8)
# text((bpctdown+0.5)[,2], ylow/100, drlabels[4:6] ,cex=0.9, pos=2, srt=90, col="black", font=2) 
text(colMeans(bpctdown) ,ylow, shortelementlabels[c(1:7,9)], cex=0.8, pos=4, srt=90, font=2, xpd=NA)
mtext("Varied TRP element(s)", side=1, line=1, cex=0.8)
arrows(bpctdown, aperm(100*rr_pctdown, c(2,1,3))[,,"0.25"], bpctdown, aperm(100*rr_pctdown, c(2,1,3))[,,"0.75"], angle=90, code=3, length=0.05)

mtext("Reduction in year 10 TB mortality with novel DR-TB regimen,\n compared to projection under current standard of care",side=3,outer=TRUE, cex=1.4, line=-1, xpd=NA) 

## DR combos
novelwide <- allnovelwide$DRDSTall_
drcombos <- numeric(0)
i <- 1; while(file.exists(paste0(location,"TRPcombos_DRDSTall_",currenttag,".",i,".csv")))
{drcombos <- rbind(drcombos, read.csv(paste0(location,"TRPcombos_DRDSTall_",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0read.csv(paste0())
outcome <- "rrdeaths"
drout <- alldrout
a <- numeric(0); for (i in 1:length(grep("^rrdeaths",colnames(drcombos)))) a[i] <- median((drcombos[,grep("^rrdeaths",colnames(drcombos))[i]]- drout$rrdeaths10 )/drout$rrdeaths10)
a[length(a)+1] <- median(((allnovelwide$DRDSTall_$rrdeaths10alloptimal - alldrout$rrdeaths10 )/alldrout$rrdeaths10))
names(a) <- c(colnames(drcombos)[grep("^rrdeaths",colnames(drcombos))], "all optimal")
par(mar=c(3,5,3,3)); 
bcols <- c("gray",topo.colors(3))
# b<-barplot(a[c(c(1,6,4,2,7,5,3)+7*rep(c(0,1,2), each=7), length(a))], las=2, ylab="mortality reduction", space = c(rep(0.1,7),1,rep(0.1,6),1,rep(0.1,6),1),
#            col=c( rep(c( bcols[1], rep(bcols[2:4], 2)), times=3), bcols[5]), angle=20, density=c( rep( c(rep(40,4), rep(100,3)), times=3), 100) ,
#            ylim=c(-0.2,0.02), xaxt='n')
b<-barplot(a[c(c(1,10,7,4,12,9,6)+12*rep(c(0,1,2), each=7))], las=2, ylab="RR TB mortality reduction", space = c(rep(0.1,7),1,rep(0.1,6),1,rep(0.1,6)),
           col=c( rep(c( bcols[1], rep(bcols[2:4], 2)), times=3)), angle=20, density=c( rep( c(rep(30,4), rep(100,3)), times=3), 100) ,
           ylim=c(-0.7,0.02), xaxt='n',  cex.lab=1.2,
           main="Impact of select characteristic combinations,\nwith and without associated increased RR detection and treatment")
text(x=b[c(4,11,18)], y=0.02, c("Current efficacy", "Intermediate efficacy", "Optimal efficacy"), xpd=NA)
legend("bottomleft", c("Baseline ~ standard of care", "Increased eligibility", "Improved adherence", "Shortened duration"),
       fill=bcols[1:4], xpd=NA, bty="n")
legend("bottomright", c("With current RR diagnosis and treatment coverage", "With universal RR diagnosis and treatment"), 
       angle=20, density = c(30,100), fill = c("black","black"), bty='n')


# Union abstract version

novelwide <- allnovelwide[["DRDSTall_"]]; drout <- droutdr <- alldrout#[1:nrow(novelwide),]

outcome <- c("rrdeaths"); y <- 10
rr_pctdown <- array(0,dim=c( length(drelementnames) , 3 , 5 )); 
dimnames(rr_pctdown) <- list("vary"=drelementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in drelementnames[-5]) 
{ rr_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, y, vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                    drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  rr_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, y,"allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                    drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  rr_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, y, vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                    drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}  

par(mar=c(6,4,2,1), mfrow=c(1,1), oma=c(0,0,2,0))

ylow <- -100
bpctdown <- barplot(height = t(100*rr_pctdown[c(1,2,6,7,3,4,8),,"0.5"]), beside = TRUE,  #width=c(1,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1),
                    ylab="% reduction (median [IQR])", 
                    main="", xlab="", las=2, cex.lab=1.1, 
                    ylim=c(1,0)*ylow, yaxt='n',
                    space=c(0,0,0,1.5,0,0,rep(c(0.5,0,0),5)),
                    col=cols, names.arg=rep("", length(drelementnames)-1))
axis(2, at=-seq(0,-ylow,by=5), labels=paste0(seq(0,-ylow,by=5),"%"),las=2,cex.axis=0.8)
text(colMeans(bpctdown)-0.5 ,ylow, c("All regimen characteristics\nvaried including\ntreatment scale-up", "Efficacy","Duration","Adherence","Barrier to resistance","Eligibility exclusions","Scale-up of\nRR TB treatment"), 
     cex=1.1, pos=4, srt=90, font=2, xpd=NA)
arrows(bpctdown, t(100*rr_pctdown[c(1,2,6,7,3,4,8),,"0.25"]), bpctdown, t(100*rr_pctdown[c(1,2,6,7,3,4,8),,"0.75"]), angle=90, code=3, length=0.05)
lpos <- legend(legend=c("minimal","intermediate","optimistic"), xjust=0.5, title="Level of varied characteristic(s) (All others held at intermediate level)", fill=cols, x=bpctdown[2,4]-0.3, y=-102,cex=1.1, xpd=NA, bty='n')
mtext("Reduction in year 10 RR-TB mortality using a novel RR-TB regimen,\n compared to projection under current standard of care",side=3,outer=TRUE, cex=1.2, line=-1, xpd=NA) 

rr_pctdown

outcome <- c("rrinc"); y <- 10
rr_pctdown <- array(0,dim=c( length(drelementnames) , 3 , 5 )); 
dimnames(rr_pctdown) <- list("vary"=drelementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in drelementnames) 
{ rr_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, y, vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                    drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  rr_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, y,"allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                    drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  rr_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, y, vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                    drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}  
rr_pctdown[,,2:4]



# RR regimen minimal and maximla baselines
novelwide <- allnovelwide[["DRDSTall_"]]; drout <- alldrout
outcome <- "rrdeaths"

drmin <- array(0,dim=c(6,5)); dimnames(drmin) <- list("element"=drelementnames[-c(1,8)], "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in drelementnames[2:7])
{
  drmin[vary,] <- quantile( 
    (allminopt_dr[ , paste0(outcome, "10allmin")] - only_dr[ , paste0(outcome, "10only",vary)] ) /
      (allminopt_dr[ , paste0(outcome, "10allmin")] - allminopt_dr[ , paste0(outcome, "10allopt")]) , c(0.025,0.25,0.5,0.75,0.975) )
}

par(mar=c(1,1,3,1), mfrow=c(1,1)) 
b <- barplot(drmin[,3], horiz = TRUE, beside=TRUE, las=2, font=2, xaxt='n',
             legend.text=shortelementlabels[-c(1,8,9)], col=rainbow(6), args.legend=list(x="topright", cex=0.8), xlim=c(-0.1,1))
mtext("Fraction of total RR TB mortality benefit of all-optimal novel RR TB regimen that is attainable by improving 
      only a single characteristic from minimal to optimal level (median [95% uncertainty range])", side=3, line=1, font=2, cex=1, xpd=NA)
text(drmin[,5]+0.05, b, paste0(round(drmin[,3]*100, 1),"%"))
arrows(drmin[,"0.025"], b, drmin[,"0.975"], b, angle=90, code=3, length=0.05)
text(0.35, b[3], "'All minimal' baseline here is a regimen with 
        ~ standard-of-care efficacy, duration, and tolerability,
     that excludes 15% due to resistance, that excludes
     10% of patients for medical reasons,
     that has low barrier to resistance,
     and that is associated with continued
     current RR TB treatment coverage levels", pos=4)

novelwide <- allnovelwide[["DRDSTall_"]]; drout <- droutdr
outcome <- "rrdeaths"

intdrmin <- array(0,dim=c(6,5)); dimnames(intdrmin) <- list("element"=drelementnames[-c(1,8)], "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in drelementnames[2:7])
{
  intdrmin[vary,] <- quantile( 
    (allminopt_dr[ , paste0(outcome, "10allmin")] - intonly_dr[ , paste0(outcome, "10only",vary)] ) /
      (allminopt_dr[ , paste0(outcome, "10allmin")] - allminopt_dr[ , paste0(outcome, "10allopt")]) , c(0.025,0.25,0.5,0.75,0.975) )
}

novelwide <- allnovelwide[["DRDSTall_"]]; drout <- droutdr
outcome <- "rrdeaths"

drmax <- array(0,dim=c(6,5)); dimnames(drmax) <- list("element"=drelementnames[-c(1,8)], "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in drelementnames[2:7])
{
  drmax[vary,] <- quantile( 
    1- ( allbut_dr[ , paste0(outcome, "10allbut",vary)] -  drout[ , paste0(outcome,"10")]) /
      (allminopt_dr[ , paste0(outcome, "10allopt")] - drout[ , paste0(outcome,"10")]) , c(0.025,0.25,0.5,0.75,0.975) )
}

par(mar=c(1,1,3,1), mfrow=c(1,1)) 
b <- barplot(drmax[,3], horiz = TRUE, beside=TRUE, las=2, font=2, xaxt='n', xlim=c(0,1),
             legend.text=shortelementlabels[-c(1,8,9)], col=rainbow(6), args.legend=list(x="topright", cex=1))#1.5, y=8, cex=0.8))
mtext("Fraction of an all-optimal novel DR regimen's mortality impact (median [95% UR])\nthat is lost when only one characteristic is reduced to its minimal value", side=3, line=1, font=2, cex=1.1, xpd=NA)
arrows(drmax[,"0.025"], b, drmax[,"0.975"], b, angle=90, code=3, length=0.05)
text(drmax[,5]+0.05, b, paste0(round(drmax[,3]*100, 1),"%"))

intdrmax <- array(0,dim=c(6,5)); dimnames(intdrmax) <- list("element"=drelementnames[-c(1,8)], "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in drelementnames[2:7])
{
  intdrmax[vary,] <- quantile( 
    1- ( intallbut_dr[ , paste0(outcome, "10allbut",vary)] -  drout[ , paste0(outcome,"10")]) /
      (allminopt_dr[ , paste0(outcome, "10allopt")] - drout[ , paste0(outcome,"10")]) , c(0.025,0.25,0.5,0.75,0.975) )
}

par(mar=c(1,1,3,1), mfrow=c(1,1)) 
b <- barplot(intdrmax[,3], horiz = TRUE, beside=TRUE, las=2, font=2, xaxt='n', xlim=c(0,1),
             legend.text=shortelementlabels[-c(1,8,9)], col=rainbow(6), args.legend=list(x="topright", cex=1))#1.5, y=8, cex=0.8))
mtext("Fraction of an all-optimal novel DR regimen's mortality impact (median [95% UR])\nthat is lost when only one characteristic is reduced to its intermediate value", side=3, line=1, font=2, cex=1.1, xpd=NA)
arrows(intdrmax[,"0.025"], b, intdrmax[,"0.975"], b, angle=90, code=3, length=0.05)
text(intdrmax[,5], b, paste0(round(intdrmax[,3]*100, 1),"%"))


alldrmin <- array(0,dim=c(12,6)); for (i in 1:6) alldrmin[((1:2)+(2*(i-1))),i] <- rbind(intdrmin[,3], drmin[,3]-intdrmin[,3])[,i]
alldrmax <- array(0,dim=c(12,6)); for (i in 1:6) alldrmax[((1:2)+(2*(i-1))),i] <- rbind(intdrmax[,3], drmax[,3]-intdrmax[,3])[,i]

par(mar=c(1,1,3,1), mfrow=c(1,1)) 
b <- barplot(alldrmax, horiz = TRUE, beside=FALSE, las=2, font=2, xaxt='n', xlim=c(0,1),
             col=rep(rainbow(6),each=2), density=c(200,25))
legend(x="topright", cex=1, legend=c(shortelementlabels[-c(1,8,9)],"Reduced to intermediate","Reduced to minimal"), fill=c(rainbow(6),"black","black"), density=c(rep(200,7),25))
mtext("Fraction of an all-optimal novel RR regimen's mortality impact (median [95% UR])\nthat is lost when only one characteristic is reduced to its intermediate or minimal value", side=3, line=1, font=2, cex=1.1, xpd=NA)
arrows(drmax[,"0.025"], b, drmax[,"0.975"], b, angle=90, code=3, length=0.05)
text(drmax[,5]+0.05, b, paste0(round(drmax[,3]*100, 1),"%"))


## FIGURE 2 ##

## 2x2 Manuscript Figure: DS/DR, min/max baselines
#  with fixed (intermediate) scale-up variables - have fixed TRPattribfrac to do this with commented if (DR/DS) lines near top
# And at baseline, assume that over the next 10 years, rif DST scaleup continues at same pace in new patients and underDST in retreatment is halved; a new regimen can make that happen faster and can also increase the final coverage.
# And ? soften "minimal" barrier to resistance parameters? - no, kept the same, will make clear minimum is a "worst case" with DST and reduces efficacy
# And include trajectories/ save all times for at least all and efficacy  - done, TRPattribfrac now saves all times
# and get intermediates for min amd opt - need to write extra script if I want to plot this
# 
# And think about how I'm going to do sensitivity analyses: PRCCs, then one-way in subset of sims? And just one-way for scaleup variables with novel regimen.)


par(mar=c(1,2,3,0), oma=c(0,4,1,1), mfrow=c(2,2)) 

# b <- barplot(dsmin[1:6,3], horiz = TRUE, beside=TRUE, las=2, font=2, xaxt='n', yaxt='n', col=rainbow(6), xlim=c(-0.1,1),
#              main="Fraction of maximal novel regimen impact that is\nachievable by a single optimized characteristic", font.main=2, xpd=NA)
b <- barplot(alldsmin[,6:1], horiz = TRUE, beside=FALSE, las=2, font=2, xaxt='n', yaxt='n', xlim=c(-0.1,0.7),
             main="Fraction of maximal novel regimen impact that is\nachievable by a single optimized characteristic", font.main=2, xpd=NA,
             col=rep(rainbow(6),each=2), density=c(25, 200))
text(dsmin[6:1,5], b, paste0(round(dsmin[6:1,3]*100, 1),"%"), pos=4)
arrows(dsmin[6:1,"0.025"], b, dsmin[6:1,"0.975"], b, angle=90, code=3, length=0.05)
text(-0.05,mean(b),"Rifampin-\nSusceptible\nregimen",pos=2, xpd=NA, cex=1.3, font=2)
text(0.3, mean(b),"A", font=2, cex=1.5)

# b <- barplot(-dsmax[1:6,3], horiz = TRUE, beside=TRUE, las=2, font=2, xaxt='n', yaxt='n', xlim=c(-1,0), col=rainbow(6),
#              main="Fraction of maximal impact lost if a\nsingle aspect of regimen is not optimized", font.main=2, xpd=NA)
b <- barplot(-alldsmax[,6:1], horiz = TRUE, beside=FALSE, las=2, font=2, xaxt='n', yaxt='n', xlim=c(-0.9,0), 
             main="Fraction of maximal impact lost if a\nsingle aspect of regimen is not optimized", font.main=2, xpd=NA,
              col=rep(rainbow(6),each=2), density=c(200,25))
arrows(-dsmax[6:1,"0.025"], b, -dsmax[6:1,"0.975"], b, angle=90, code=3, length=0.05)
text(-dsmax[6:1,5], b, paste0(round(dsmax[6:1,3]*100, 1),"%"),xpd=NA, pos=2)
# text(-1.25,b+0.4,
#   c("Efficacy (94-99% durable cure*)","Barrier to resistance (0-5% acquire new resistance)","Preexisting novel-regimen resistance (prevalence 0-10%)","Medical Contraindications (0-11% excluded)", "Duration (2-6 months)", "Tolerability (adherence improves 0-50%)"),pos=1, xpd=NA)
text(-0.4, mean(b),"B", font=2, cex=1.5)

# b <- barplot(drmin[1:6,3], horiz = TRUE, beside=TRUE, las=2, font=2, xaxt='n', yaxt='n', col=rainbow(6), xlim=c(-0.1,1))
b <- barplot(alldrmin[,6:1], horiz = TRUE, beside=FALSE, las=2, font=2, xaxt='n', yaxt='n', xlim=c(-0.1,0.7),
             col=rep(rainbow(6),each=2), density=c(25, 200))
text(drmin[6:1,5], b, paste0(round(drmin[6:1,3]*100, 1),"%"), pos=4)
arrows(drmin[6:1,"0.025"], b, drmin[6:1,"0.975"], b, angle=90, code=3, length=0.05)
text(-0.05,mean(b),"Rifampin-\nResistant\nregimen",pos=2, xpd=NA, cex=1.3, font=2)
text(0.3, mean(b),"C", font=2, cex=1.5)

# b <- barplot(-drmax[1:6,3], horiz = TRUE, beside=TRUE, las=2, font=2, xaxt='n', yaxt='n', xlim=c(-1,0), col=rainbow(6))
b <- barplot(-alldrmax[,6:1], horiz = TRUE, beside=FALSE, las=2, font=2, xaxt='n', yaxt='n', xlim=c(-0.9,0),
             col=rep(rainbow(6),each=2), density=c(200,25))
arrows(-drmax[6:1,"0.025"], b, -drmax[6:1,"0.975"], b, angle=90, code=3, length=0.05)
text(-drmax[6:1,5], b, paste0(round(drmax[6:1,3]*100, 1),"%"), pos=2)
# text(-1.25,b+0.4,
#      c("Efficacy (76-94% durable cure*)","Barrier to resistance (0.8-10% acquire new resistance)","Preexisting resistance (prevalence 0-15%)","Medical Contraindications (0-11% excluded)", "Duration (6-20 months)", "Tolerability (adherence improves 0-50%)"),pos=1, xpd=NA)
text(-0.4, mean(b),"D", font=2, cex=1.5)

legend(x=-1.1,y=mean(b)*3.5,,xpd=NA, cex=1, legend=c(shortelementlabels[-c(1,8,9)],"Minimal vs. intermediate value","Intermediate vs. optimal value"), fill=c(rainbow(6),"black","black"), density=c(rep(200,6),25, 200))


## FIGURE 3 ##

###########
# Trajectories figure
# there's no "median trajectory" so I'll have to take the median for each intervention and each time point separately
# need to use a min and opt with int scaleup, so will redo to save these at all time points and reload on 4/2
# and if using RR, need to update axis labels and manuscript text

novelwide <- novelwide3
outcome <- "rrdeaths"
times <- seq(0,10,by=1); ls <- c("allmin", "onlyeff", "allbuteff", "allopt")
traj <- array(0, dim=c(length(times),4,5)); dimnames(traj) <- list("t"=times, "level"=ls, "q"=c(0.025,0.25,0.5,0.75,0.975))
for (t in times) 
{
  traj[t+1,"allmin",] <- quantile(allminopt_dr[,colnames(allminopt_dr)==paste0(outcome, t, "allmin")], c(0.025,0.25,0.5,0.75,0.975))
  traj[t+1,"allopt",] <- quantile(allminopt_dr[,colnames(allminopt_dr)==paste0(outcome, t, "allopt")], c(0.025,0.25,0.5,0.75,0.975))
  traj[t+1,"onlyeff",] <- quantile(only_dr[,colnames(only_dr)==paste0(outcome, t, "onlyefficacy")], c(0.025,0.25,0.5,0.75,0.975))
  traj[t+1,"allbuteff",] <- quantile(allbut_dr[,colnames(allbut_dr)==paste0(outcome, t, "allbutefficacy")], c(0.025,0.25,0.5,0.75,0.975))
#   traj[t+1,"onlyeff",] <- quantile(only_dr[,colnames(only_dr)==paste0(outcome, t, "onlyefficacy")], c(0.025,0.25,0.5,0.75,0.975))
#   traj[t+1,"allbuteff",] <- quantile(allbut_dr[,colnames(allbut_dr)==paste0(outcome, t, "allbutefficacy")], c(0.025,0.25,0.5,0.75,0.975))
}
par(mar=c(4,4,1,1), oma=c(0,0,0,0), mfrow=c(1,1))
plot(0:10, traj[,"allmin","0.5"], ylim=c(0,1.4), type='l', xlab="Years after novel RR TB regimen's introduction", ylab="Annual RR TB mortality, median projection",
# plot(0:10, traj[,"allmin","0.5"], ylim=c(0,40), type='l', xlab="Years after novel RS TB regimen's introduction", ylab="Annual TB mortality, median projection",
          lwd=2, col="red",lty=1, bty='l')
points(0:10, traj[,"onlyeff","0.5"], type='l', col='red',lwd=2, lty=2)
points(0:10, traj[,"allbuteff","0.5"], type='l', col='green',lwd=2, lty=2)
points(0:10, traj[,"allopt","0.5"], type='l', col='green',lwd=2, lty=1)
# segments(x0=10, y0=20, y1=35, lty=3, col="black", lwd=1); text(10,18,"Primary analyses\nperformed 10 years\nafter novel regimen's\nintroduction", pos=2)
segments(x0=10, y0=0.5, y1=1.5, lty=3, col="black", lwd=1); text(10,0.5,"Primary analyses\nperformed 10 years\nafter novel regimen's\nintroduction", pos=2)

legend(x=0.2,y=0.01,xjust=0,yjust=0, legend= c("RR TB regimen minimal for all characteristics", "Improved efficacy only", "All characteristics improved except efficacy", "All characteristics improved"),
       col=c("red","red","green","green"), lty=c(1,2,2,1),lwd=2, cex=1)

arrows(9, traj["9","allbuteff","0.5"], 9, traj["9","allopt","0.5"], angle=45, code=3, length=0.1); 
          text (9,mean(traj["9",c("allopt","onlyeff"),"0.5"]),"Figure 2C", pos=2)
arrows(9.5, mean(traj[c("10","9"),"onlyeff","0.5"]), 9.5, mean(traj[c("10","9"),"allmin","0.5"]), angle=45, code=3, length=0.1); 
          text (9.5,mean(traj[c("10","9"),c("onlyeff","allmin"),"0.5"]),"Figure 2D", pos=2)
# add curved arrows after generating figure?




############ duration
# DS regimen:
# resource use: barplots as above but for outcomes of diagnoses, DSTs (rif and novel in same plot), and rxmonths (all 3 in same plot)
# switch to al-min baseline
novelwide <- novelwide1 <- allnovelwide[["DSDSTall_rDSTall."]]; drout <- alldrDST

cumresource <- array(0,dim=c( length(elementnames[2:7]) , 3 , 3)); resource <- array(0,dim=c( length(elementnames[2:7]) , 3 , 3,5)); 
dimnames(cumresource) <- list("vary"=elementnames[2:7], "level"=c("minimal", "intermediate", "optimal"), "reg"=c("First-line","Second-line","Novel")); dimnames(resource) <- list("vary"=elementnames[2:7], "level"=c("minimal", "intermediate", "optimal"), "reg"=c("First-line","Second-line","Novel"), "q"=c(0.025,0.25,0.5,0.75,0.975))

for (vary in dselementnames[2:7]) for (nreg in 1:3) 
{ 
  outcome <- c("rxtime_s","rxtime_r","rxtime_n", "dxs")[nreg]
  resource[vary,1,nreg,] <- quantile(allminopt_ds[ , paste0(outcome, "10allmin")],c(0.025,0.25,0.5,0.75,0.975))
  resource[vary,2,nreg,] <- quantile(intonly_ds[ , paste0(outcome, "10only",vary)],c(0.025,0.25,0.5,0.75,0.975))
  resource[vary,3,nreg,] <- quantile(only_ds[ , paste0(outcome, "10only",vary)],c(0.025,0.25,0.5,0.75,0.975))
}  #   resource[vary,1,nreg,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"minimal")],c(0.025,0.25,0.5,0.75,0.975))
#   resource[vary,2,nreg,] <- quantile(novelwide[ , paste0(outcome, "10allintermediate")],c(0.025,0.25,0.5,0.75,0.975))
#   resource[vary,3,nreg,] <- quantile(novelwide[ , paste0(outcome, "10", vary,"optimal")],c(0.025,0.25,0.5,0.75,0.975))
#   for (t in 1:10)
#   {
#     cumresource[vary,1,nreg] <- cumresource[vary,1,nreg] + median(novelwide[ , paste0(outcome, t, vary,"minimal")])
#     cumresource[vary,2,nreg] <- cumresource[vary,2,nreg] + median(novelwide[ , paste0(outcome, t,"allintermediate")])
#     cumresource[vary,3,nreg] <- cumresource[vary,3,nreg] + median(novelwide[ , paste0(outcome, t, vary,"optimal")])
#   }
# }  

## FIGURE 4 ##

par(mar=c(3,6,4,1), mfcol=c(1,1), oma=c(0,0,0,0)) 
bres <- barplot(12*cbind(aperm(resource,c(1,3,2,4))["efficacy",c(2,1,3),c(1,3),3], aperm(resource,c(1,3,2,4))["exclusions",c(2,1,3),3,3],aperm(resource,c(1,3,2,4))["duration",c(2,1,3),3,3]), beside = FALSE, 
                cex.lab=1.2, ylim=c(0,600), #main="Treatment provided", 
                col=rev(blues), ylab="Patient-months of treatment\nper 100,000 population in year 10, by regimen", xaxt='n')
text(bres, -50, c("All-minimal\nRS regimen","Optimized\nEfficacy", 
                             "Optimized\nEligibility","Optimized\nDuration"), font=2, xpd=NA)
legend(x = bres[3]+0.5, y=max(12*rowSums(resource["duration",,1:3,3]))+30, xjust=0.5, yjust=0, fill=blues,
       c("Novel DS regimen","MDR regimen","Standard DS regimen")[c(1,3,2)], xpd=NA)
# text((bres[c(1,3,5)]+bres[c(2,4,6)])/2, -6*14, c("Efficacy", "Contraindications","Duration"), cex=1.1, pos=1, xpd=NA)


# diagnostic and other resources

tests <- array(0,dim=c( length(dselementnames) , 3 , 3, 5 ));
dimnames(tests) <- list("vary"=dselementnames, "level"=c("minimal", "intermediate", "optimal"), "test"=c("Novel regimen DSTs performed","Rifampin DSTs performed","Treatment courses initiated"),"q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in dselementnames)  for (ntest in 1:3)
{ test <- c("nDSTs","rDSTs","dxs")[ntest]
  tests[vary,1,ntest,] <- quantile((1-novelwide$initialloss_s)*novelwide[ , paste0(test, "10", vary,"minimal")], c(0.025,0.25,0.5,0.75,0.975))
  tests[vary,2,ntest,] <- quantile((1-novelwide$initialloss_s)*novelwide[ , paste0(test, "10allintermediate")], c(0.025,0.25,0.5,0.75,0.975))
  tests[vary,3,ntest,] <- quantile((1-novelwide$initialloss_s)*novelwide[ , paste0(test, "10", vary,"optimal")], c(0.025,0.25,0.5,0.75,0.975))
}  

par(mar=c(5,6,2,1), mfcol=c(1,1), oma=c(0,0,0,0)) 
bup <- barplot(c(aperm(tests,c(1,3,2,4))["efficacy",1,,3], aperm(tests,c(1,3,2,4))["companion",1,,3],aperm(tests,c(1,3,2,4))["duration",1,,3]), beside = TRUE, 
               space=c(0.75,0.25,0.25), ylim=c(0,70), cex.lab=1.2, main="Diagnostic testing",
               col=cols, ylab="Novel regimen DSTs\nperformed, year 10", xlab="")
text(bup[c(2,5,8)], -10, c("% Durably Cured", "Baseline novel-\nregimen resistance","Duration"), cex=1.1, pos=1, xpd=NA)
arrows(bup, aperm(tests, c(2,1,3,4))[,2:4,1,"0.25"], bup, aperm(tests, c(2,1,3,4))[,2:4,1,"0.75"], angle=90, code=3, length=0.05)


# DR regimen:
novelwide <- allnovelwide[["DRDSTall_"]]; drout <- alldrout[1:nrow(novelwide),]

cumresource <- array(0,dim=c( length(drelementnames) , 3 , 3)); resource <- array(0,dim=c( length(drelementnames) , 3 , 3,5)); 
dimnames(cumresource) <- list("vary"=drelementnames, "level"=c("minimal", "intermediate", "optimal"), "reg"=c("Standard DS regimen","Standard DR regimen","Novel regimen")); dimnames(resource) <- list("vary"=drelementnames, "level"=c("minimal", "intermediate", "optimal"), "reg"=c("First-line","Second-line","Novel DR"), "q"=c(0.025,0.25,0.5,0.75,0.975))

for (vary in drelementnames) for (nreg in 1:3) 
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

par(mar=c(7,6,4,1), mfcol=c(1,1), oma=c(0,0,0,0)) 
bres <- barplot(12*cbind(aperm(resource,c(1,3,2,4))["efficacy",2:3,,3], aperm(resource,c(1,3,2,4))["companion",2:3,,3],aperm(resource,c(1,3,2,4))["duration",2:3,,3]), beside = FALSE, 
                space=c(0.75,0.25,0.25), cex.lab=1.2, main="DR-specific treatment provided", ylim=c(0,40),
                col=rev(blues), ylab="Patient-months of treatment\nin year 10, by regimen")
legend(x = bres[3], y=max(12*rowSums(resource["duration",,2:3,3]))-4, xjust=0.5, yjust=0, fill=blues[2:3],
       c("Novel DR regimen","Standard DR regimen"), xpd=NA)
text(bres[7:9], -4, c("20mo","9mo","6mo"), cex=0.9, pos=1, xpd=NA)
text(bres[c(2,5,8)], -6, c("% Durably Cured", "Baseline novel-\nregimen resistance","Duration"), cex=1.1, pos=1, xpd=NA)

tests <- array(0,dim=c( length(drelementnames) , 3 , 3, 5 ));
dimnames(tests) <- list("vary"=drelementnames, "level"=c("minimal", "intermediate", "optimal"), "test"=c("Novel regimen DSTs performed","Rifampin DSTs performed","Treatment courses initiated"),"q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in drelementnames)  for (ntest in 1:3)
{ test <- c("nDSTs","rDSTs","dxs")[ntest]
  tests[vary,1,ntest,] <- quantile((1-novelwide$initialloss_s)*novelwide[ , paste0(test, "10", vary,"minimal")], c(0.025,0.25,0.5,0.75,0.975))
  tests[vary,2,ntest,] <- quantile((1-novelwide$initialloss_s)*novelwide[ , paste0(test, "10allintermediate")], c(0.025,0.25,0.5,0.75,0.975))
  tests[vary,3,ntest,] <- quantile((1-novelwide$initialloss_s)*novelwide[ , paste0(test, "10", vary,"optimal")], c(0.025,0.25,0.5,0.75,0.975))
}  

par(mar=c(5,6,2,1), mfcol=c(1,1), oma=c(0,0,0,0)) 
bup <- barplot(c(aperm(tests,c(1,3,2,4))["efficacy",1,,3], aperm(tests,c(1,3,2,4))["companion",1,,3],aperm(tests,c(1,3,2,4))["duration",1,,3]), beside = TRUE, 
               space=c(0.75,0.25,0.25), ylim=c(0,1), cex.lab=1.2, main="Up-front costs:\nDiagnostic testing for novel DR regimen",
               col=cols, ylab="Novel DR regimen DSTs performed\n(year 10, per 100K population)", xlab="")
text(bup[c(2,5,8)], -0.2, c("% Durably Cured", "Baseline novel-\nregimen resistance","Duration"), cex=1.1, pos=1, xpd=NA)
arrows(bup, aperm(tests, c(2,1,3,4))[,2:4,1,"0.25"], bup, aperm(tests, c(2,1,3,4))[,2:4,1,"0.75"], angle=90, code=3, length=0.05)


##############
# Sensitivity analysis, total regimen impact
library("sensitivity")
prccall <- pcc(X = novelwide1[,41:74], y= ((alldrDST[,"tbdeaths10"] - novelwide1[,"tbdeaths10allintermediate"])/ alldrDST[,"tbdeaths10"] ), rank=TRUE ) 
prccallres <- pcc(X = novelwide3[,41:74], y= ((alldrout[,"rrdeaths10"] - novelwide3[,"rrdeaths10allintermediate"]) /alldrout[,"rrdeaths10"] ), rank=TRUE ) 
cbind(rownames(prccall$PRCC), prccall$PRCC)[rev(order(abs(prccall$PRCC))),]

## SUPPLEMENTAL FIGURE S[all]##
o <- rev(rev(order(abs(prccall$PRCC$original)))[1:10])
ores <- rev(rev(order(abs(prccallres$PRCC$original)))[1:10])
display <- c(ores[!(ores %in% o)],o)

par(mfrow=c(1,2),mar=c(4,2,2,1), oma=c(0,25,3,0))
b <- barplot(prccall$PRCC$original[display], horiz=TRUE, beside=TRUE, space=c(0.5,0),
             names.arg=longparamnames[display], 
             axes=FALSE,axisnames=FALSE, xlim=c(-1,1),
             las=1, cex.names=0.8, cex.axis=0.8,cex.lab=0.8, main="Novel DS regimen", cex.main=1 )
axis(1, at=0.2*c(-4:4),cex.axis=0.8,cex.lab=0.9)
axis(2, labels=longparamnames[display],at=b, las=1, cex.axis=0.8, cex.lab=2, outer = TRUE, tick = FALSE, xpd=NA)
mtext("Parameter", side=2,line=24, outer=TRUE)
b <- barplot(prccallres$PRCC$original[display], horiz=TRUE, beside=TRUE, space=c(0.5,0),
             names.arg=longparamnames[display], 
             axes=FALSE,axisnames=FALSE, xlim=c(-1,1),
             las=1, cex.names=0.8, cex.axis=0.8,cex.lab=0.8, main="Novel DR regimen", cex.main=1 )
axis(1, at=0.2*c(-4:4),cex.axis=0.8,cex.lab=0.9)
mtext("Partial rank correlation of model parameters with the\nmortality impact of an all-intermediate novel regimen", side=3,cex=1.2,outer=TRUE)



# # sensitivity analysis: ltfu at 2 months (haven't run for 20160313p)
# ltfu2mo <- read.csv("TRPwideoutput_DSDSTall_ltfu2mo.rDSTall.India_20160201.1.csv", header=TRUE)
# drltfu <- screendrout("DRcalibration_ltfu2mo.rDSTall.India_20160201.1.csv", tolerance=1.5)
# outcome <- "tbdeaths"
# down <- array(0,dim=c( length(elementnames) , 3 , 5 )); 
# dimnames(down) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
# for (vary in elementnames)
# { down[vary,1,] <- quantile((ltfu2mo[ , paste0(outcome, "10", vary,"minimal")] - drltfu[ , paste0(outcome,"10")] )/
#                               drltfu[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
#   down[vary,2,] <- quantile((ltfu2mo[ , paste0(outcome, "10allintermediate")] - drltfu[ , paste0(outcome,"10")] )/
#                               drltfu[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
#   down[vary,3,] <- quantile((ltfu2mo[ , paste0(outcome, "10", vary,"optimal")] - drltfu[ , paste0(outcome,"10")] )/
#                               drltfu[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
# }  
# 
# par(mar=c(2,4,3,1), mfrow=c(1,1), oma=c(0,0,0,0))
# bpctdown <- barplot(height = 100*aperm(down, c(2,1,3))[,,"0.5"], beside = TRUE, 
#                     ylab="% reduction (median [IQR])",
#                     main="With all losses to follow up occurring at 2 months\n(to maximize the impact of regimen duration)",
#                     xlab="", las=2, cex.lab=1, 
#                     ylim=c(150*final_pctdown[1,3,"0.5"], 2), yaxt='n',
#                     #                     legend=c("minimal","intermediate","optimal"),
#                     #                     args.legend=list(title="Level of varied element(s)", x=25,y=-15,cex=0.9),
#                     space=c(0,0,0,1.5,0,0,rep(c(0.5,0,0),5)),
#                     col=cols, names.arg=rep("", length(elementnames)))
# axis(2, at=seq(-25,0,by=5), labels=paste0(seq(-25,0,by=5),"%"),las=2,cex.axis=0.8)
# # text((bpctdown+0.4)[,2], -0.5, dslabels[4:6] ,cex=0.9, pos=2, srt=90, col="black", font=2) 
# text(colMeans(bpctdown) ,-23, elementlabels, cex=0.9, pos=1, srt=90, font=1, xpd=NA)
# mtext("Varied TRP element(s)", side=1, line=1, cex=0.9)
# arrows(bpctdown, aperm(100*down, c(2,1,3))[,,"0.25"], bpctdown, aperm(100*down, c(2,1,3))[,,"0.75"], angle=90, code=3, length=0.05, xpd=NA)
# # arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)


# sensitivity analysis: PRCCs for each parameter
prcc$duration <- pcc(X = novelwide1[,41:74], y= (allminopt[,"tbdeaths10durationminimal"] - novelwide1[,"tbdeaths10durationoptimal"] )/ (novelwide1[,"tbdeaths10allminimal"] - novelwide1[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$durationres <- pcc(X = novelwide3[,41:74], y= (novelwide3[,"tbdeaths10durationminimal"] - novelwide3[,"tbdeaths10durationoptimal"] )/ (novelwide3[,"tbdeaths10allminimal"] - novelwide3[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$duration <- pcc(X = novelwide1[,41:74], y= (allminopt[,"tbdeaths10durationminimal"] - novelwide1[,"tbdeaths10durationoptimal"] )/ (novelwide1[,"tbdeaths10allminimal"] - novelwide1[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$durationres <- pcc(X = novelwide3[,41:74], y= (novelwide3[,"tbdeaths10durationminimal"] - novelwide3[,"tbdeaths10durationoptimal"] )/ (novelwide3[,"tbdeaths10allminimal"] - novelwide3[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$duration <- pcc(X = novelwide1[,41:74], y= (allminopt[,"tbdeaths10durationminimal"] - novelwide1[,"tbdeaths10durationoptimal"] )/ (novelwide1[,"tbdeaths10allminimal"] - novelwide1[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$durationres <- pcc(X = novelwide3[,41:74], y= (novelwide3[,"tbdeaths10durationminimal"] - novelwide3[,"tbdeaths10durationoptimal"] )/ (novelwide3[,"tbdeaths10allminimal"] - novelwide3[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$duration <- pcc(X = novelwide1[,41:74], y= (allminopt[,"tbdeaths10durationminimal"] - novelwide1[,"tbdeaths10durationoptimal"] )/ (novelwide1[,"tbdeaths10allminimal"] - novelwide1[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$durationres <- pcc(X = novelwide3[,41:74], y= (novelwide3[,"tbdeaths10durationminimal"] - novelwide3[,"tbdeaths10durationoptimal"] )/ (novelwide3[,"tbdeaths10allminimal"] - novelwide3[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$duration <- pcc(X = novelwide1[,41:74], y= (allminopt[,"tbdeaths10durationminimal"] - novelwide1[,"tbdeaths10durationoptimal"] )/ (novelwide1[,"tbdeaths10allminimal"] - novelwide1[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$durationres <- pcc(X = novelwide3[,41:74], y= (novelwide3[,"tbdeaths10durationminimal"] - novelwide3[,"tbdeaths10durationoptimal"] )/ (novelwide3[,"tbdeaths10allminimal"] - novelwide3[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$duration <- pcc(X = novelwide1[,41:74], y= (allminopt[,"tbdeaths10durationminimal"] - novelwide1[,"tbdeaths10durationoptimal"] )/ (novelwide1[,"tbdeaths10allminimal"] - novelwide1[,"tbdeaths10alloptimal"] ), rank=TRUE)
prcc$durationres <- pcc(X = novelwide3[,41:74], y= (novelwide3[,"tbdeaths10durationminimal"] - novelwide3[,"tbdeaths10durationoptimal"] )/ (novelwide3[,"tbdeaths10allminimal"] - novelwide3[,"tbdeaths10alloptimal"] ), rank=TRUE)

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

##############
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

bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,"companion","0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality (median [IQR])", 
                    xlab="Baseline prevalence of\ncompanion-drug resistance", names.arg=c("10%","3%","0%"), cex.lab=1, main="With DST for\nnovel DS regimen",
                    ylim=c(-15,2), col=cols)
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,"companion","0.25"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,"companion","0.75"], angle=90, code=3, length=0.05)


bpctdown <- barplot(height = 100*aperm(final_pctdown2, c(2,1,3))[,"companion","0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality (median [IQR])", 
                    xlab="Baseline prevalence of\ncompanion-drug resistance", names.arg=c("10%","3%","0%"), cex.lab=1, main="Without DST for\nnovel DS regimen",
                    ylim=c(-15,2), col=cols)
arrows(bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,"companion","0.25"], bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,"companion","0.75"], angle=90, code=3, length=0.05, xpd=NA)



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
arrows(bpctdown, aperm(100*r_pctdown2, c(2,1,3))[,4,"0.25"], bpctdown, aperm(100*r_pctdown2, c(2,1,3))[,4,"0.75"], angle=90, code=3, length=0.05, xpd=NA)


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
o2 <- rev(rev(order(abs(prcc$companionnodst$PRCC$original)))[1:10])
ores2 <- rev(rev(order(abs(prcc$companionnodst_dr$PRCC$original)))[1:10])
display <- c(o2[!o2 %in% c(ores, o)], ores2[!ores2 %in% c(o2,ores,o)], ores[!(ores %in% o)],o)

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

bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,"barrier","0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality (median [IQR])", 
                    xlab="Per-treatment probability\nof acquired novel resistance", names.arg=c("1/20","1/125","0"), cex.lab=1, main="With DST for\nnovel DS regimen",
                    ylim=c(-15,2), col=cols)
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,"barrier","0.25"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,5,"0.75"], angle=90, code=3, length=0.05)


bpctdown <- barplot(height = 100*aperm(final_pctdown2, c(2,1,3))[,"barrier","0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality (median [IQR])", 
                    xlab="Per-treatment probability\nof acquired novel resistance", names.arg=c("1/20","1/125","0"), cex.lab=1, main="Without DST for\nnovel DS regimen",
                    ylim=c(-15,2), col=cols)
arrows(bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,5,"0.25"], bpctdown, aperm(100*final_pctdown2, c(2,1,3))[,"barrier","0.75"], angle=90, code=3, length=0.05)


#DR:
par(mar=c(5,5,3,1), mfrow=c(1,2), oma=c(0,0,0,0))

bpctdown <- barplot(height = 100*aperm(rr_pctdown, c(2,1,3))[,"barrier","0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 rifampin-resistant\nTB mortality (median [IQR])", 
                    xlab="Per-treatment probability\nof acquired novel resistance", names.arg=c("1/10","1/20","1/125"), cex.lab=1, 
                    main="With DST for\nnovel DR regimen",
                    ylim=c(-40,2), col=cols)
arrows(bpctdown, aperm(100*rr_pctdown, c(2,1,3))[,"barrier","0.25"], bpctdown, aperm(100*rr_pctdown, c(2,1,3))[,"barrier","0.75"], angle=90, code=3, length=0.05)

bpctdown <- barplot(height = 100*aperm(r_pctdown2, c(2,1,3))[,"barrier","0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 rifampin-resistant\nTB mortality(median [IQR])", 
                    xlab="Per-treatment probability\nof acquired novel resistance", names.arg=c("1/10","1/20","1/125"), cex.lab=1, 
                    main="Without DST for\nnovel DR regimen",
                    ylim=c(-40,2), col=cols)
arrows(bpctdown, aperm(100*r_pctdown2, c(2,1,3))[,"barrier","0.25"], bpctdown, aperm(100*r_pctdown2, c(2,1,3))[,"barrier","0.75"], angle=90, code=3, length=0.05)


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
o2 <- rev(rev(order(abs(prcc$barriernodst$PRCC$original)))[1:10])
ores2 <- rev(rev(order(abs(prcc$barriernodst_dr$PRCC$original)))[1:10])
display <- c(o2[!o2 %in% c(ores, o)], ores2[!ores2 %in% c(o2,ores,o)], ores[!(ores %in% o)],o)

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

rtallynames <- readRDS("rtallynames_20160111.RDS")

RtrajDST <- read.csv("Resistance_DSDSTall_rDSTall.India_20160201.idr1.csv")
for (i in 2:5) RtrajDST <- rbind(RtrajDST, read.csv(paste0("Resistance_DSDSTall_rDSTall.India_20160201.idr",i,".csv")))
drout <- rbind(alldrDST[alldrDST$idr == 1,], alldrDST[alldrDST$idr == 2,], alldrDST[alldrDST$idr == 3,],
               alldrDST[alldrDST$idr == 4,], alldrDST[alldrDST$idr == 5,])

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
 if (b==3) mtext(paste0(c("Unfavorable (10%)","Intermediate (3%)","No")[c]," baseline\ncompanion-drug resistance"),side=3, line =1, cex=1,font=2)
 if (c==3) mtext(paste0(c("Minimal","Intermediate","Optimal")[b]," barrier to resistance\n"),side=2, cex=1,line=4, font=2)
                        #,c("(1/20 treatments","(1/125 treatments", "(no acquired resistance")[b],"\nif iniitally regimen- susceptible)"
#   mtext(paste0(c("Unfavorable (10%)","Optimal (3%)","No")[c]," baseline companion-resistance,\nand ",levels[b]," barrier to resistance\n(",c("1/20 treatments)","1/125 treatments)", "no acquired resistance if iniitally fully susceptible")[b]),
#         side=3,line=1,xpd=NA, cex=0.7)
  mtext("Year",side=1,line=2,cex=0.7)
  axis(side=2,at = c(seq(0,0.3,by=0.05)), labels = paste0(seq(0,30,by=5),"%"),las=2)
points(0:10, rtraj[b,c,"nprev",3,], type='l', lwd=2,col="purple")
points(0:10, rtraj[b,c,"cprev",3,], type='l', lwd=2,col="green")
}
mtext("With DST for novel regimen",side=3,outer=TRUE,line=3,cex=1.3)

Rtraj <- read.csv("Resistance_DSDSTnone_rDSTall.India_20160201.idr1.csv")
for (i in 2:5) Rtraj <- rbind(Rtraj, read.csv(paste0("Resistance_DSDSTnone_rDSTall.India_20160201.idr",i,".csv")))
drout <- rbind(alldrDST[alldrDST$idr == 1,], alldrDST[alldrDST$idr == 2,], alldrDST[alldrDST$idr == 3,],
               alldrDST[alldrDST$idr == 4,], alldrDST[alldrDST$idr == 5,])

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
  if (b==3) mtext(paste0(c("Unfavorable (10%)","Intermediate (3%)","No")[c]," baseline\ncompanion-drug resistance"),side=3, line =1, cex=1,font=2)
 if (c==3) mtext(paste0(c("Minimal","Intermediate","Optimal")[b]," barrier to resistance\n"),side=2, cex=1,line=3, font=2)
                        #c("(1 in 20 acquire resistance","(1 in 125 acquire resistance", "(no acquired resistance")[b],"\nif iniitally regimen- susceptible)"),side=2, cex=0.8,line=4, font=2) #   mtext(paste0(c("Unfavorable (10%)","Optimal (3%)","No")[c]," baseline companion-resistance,\nand ",levels[b]," barrier to resistance\n(",c("1/20 treatments)","1/125 treatments)", "no acquired resistance if iniitally fully susceptible")[b]),
         #side=3,line=1,xpd=NA, cex=1)
 mtext("Year",side=1,line=2,cex=0.7)
 axis(side=2,at = c(seq(0,0.3,by=0.05)), labels = paste0(seq(0,30,by=5),"%"),las=2)
 points(0:10, rtraj[b,c,"nprev",3,], type='l', lwd=2,col="purple")
 points(0:10, rtraj[b,c,"cprev",3,], type='l', lwd=2,col="green")
}
mtext("Without DST for novel regimen",side=3,outer=TRUE,line=3,cex=1.3)




### varying both companion and barrier

RtrajDST <- read.csv("Resistance_DSDSTall_rDSTall.India_20160201.idr1.csv")
for (i in 2:5) RtrajDST <- rbind(RtrajDST, read.csv(paste0("Resistance_DSDSTall_rDSTall.India_20160201.idr",i,".csv")))
drout <- rbind(alldrDST[alldrDST$idr == 1,], alldrDST[alldrDST$idr == 2,], alldrDST[alldrDST$idr == 3,],
               alldrDST[alldrDST$idr == 4,], alldrDST[alldrDST$idr == 5,])

resDST <- array(0,dim=c( 3 , 3 , 5 ));
dimnames(resDST) <- list("companion"=levels, "barrier"=levels, "q"=c(0.025,0.25,0.5,0.75,0.975))
outcome <- "tbdeaths"
for (c in levels) for (b in levels)
{ 
  resDST[c,b,] <- quantile((RtrajDST[ , paste0(outcome, "10companion", c,"barrier",b)] - drout[ , paste0(outcome,"10")] )/
                                        drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}  


Rtraj <- read.csv("Resistance_DSDSTnone_rDSTall.India_20160201.idr1.csv")
for (i in 2:5) Rtraj <- rbind(Rtraj, read.csv(paste0("Resistance_DSDSTnone_rDSTall.India_20160201.idr",i,".csv")))
drout <- rbind(alldrDST[alldrDST$idr == 1,], alldrDST[alldrDST$idr == 2,], alldrDST[alldrDST$idr == 3,],
               alldrDST[alldrDST$idr == 4,], alldrDST[alldrDST$idr == 5,])
resnoDST <- array(0,dim=c( 3 , 3 , 5 ));
dimnames(resnoDST) <- list("companion"=levels, "barrier"=levels, "q"=c(0.025,0.25,0.5,0.75,0.975))
outcome <- "tbdeaths"
for (c in levels) for (b in levels)
{ 
  resnoDST[c,b,] <- quantile((Rtraj[ , paste0(outcome, "10companion", c,"barrier",b)] - drout[ , paste0(outcome,"10")] )/
                             drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}  


par(mar=c(5,5,3,1), mfrow=c(1,2), oma=c(0,0,0,0))

bpctdown <- barplot(height = 100*diag(resDST[,,"0.5"]), beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality (median [IQR])", 
                    xlab="Varying both companion resistance prevalence\nand barrier to acquired resistance",
                    names.arg=levels, cex.lab=1, main="With DST for\nnovel DS regimen",
                    ylim=c(-15,10), col=cols)
arrows(bpctdown, 100*diag(resDST[,,"0.25"]), bpctdown, 100*diag(resDST[,,"0.75"]), angle=90, code=3, length=0.05)

bpctdown <- barplot(height = 100*diag(resnoDST[,,"0.5"]), beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality (median [IQR])", 
                    xlab="Varying both companion resistance prevalence\nand barrier to acquired resistance",
                    names.arg=levels, cex.lab=1, main="Without DST for\nnovel DS regimen",
                    ylim=c(-15,10), col=cols)
arrows(bpctdown, 100*diag(resnoDST[,,"0.25"]), bpctdown, 100*diag(resnoDST[,,"0.75"]), angle=90, code=3, length=0.05, xpd=NA)


# Projections without novel regimen #

rrfracs <- apply(alldrout[,paste0("rrfrac",c("","10"))], 2, quantile, c(0.025,0.25,0.5,0.75,0.975))
colnames(rrfracs) <- c("year 0", "year 10 without novel regimen")
rrfracs


rrfracsDST <- apply(alldrDST[,paste0("rrfrac",c("","10"))], 2, quantile, c(0.025,0.25,0.5,0.75,0.975))
colnames(rrfracsDST) <- c("year 0", "year 10 without novel regimen")
rrfracsDST

############################
# exclusions

# India exclusions

iexc <- rbind(read.csv("Exclusions_DRDSTall_India_20160201.idr1.csv"), read.csv("Exclusions_DRDSTall_India_20160201.idr2.csv"))
# colnames(iexc) <- colnames(read.csv("Exclusions_DSDSTall_rDSTall.India_20160111.idr1.csv"))
drout <- rbind(alldrout[alldrout$idr == 1,], alldrout[alldrout$idr == 2,])
novelnow <- allnovelwide[["DRDSTall_"]]; novelnow <- rbind(novelnow[novelnow$idr==1,], novelnow[novelnow$idr==2,])
outcome <- "rrdeaths"
idown <- array(0,dim=c(  3, 3, 2, 5 ));
dimnames(idown) <- list( "efflevel"=levels, "exclevel"=levels, "HIV"=c("HIV","nonHIV"), "q"=c(0.025,0.25,0.5,0.75,0.975))

for (efficacy in levels) for (exclusions in levels) for (H in c("HIV","nonHIV"))
{  idown[efficacy,exclusions, H,] <- quantile((iexc[ , paste0(outcome, "10", H,"exclusions",exclusions,"efficacy",efficacy)] - novelnow[ , paste0(outcome,"10allminimal")] )/
                               novelnow[ , paste0(outcome,"10allminimal")], c(0.025,0.25,0.5,0.75,0.975))
}

par(mar=c(5,5,1,1), mfrow=c(1,2), oma=c(0,0,2,0))

inoH <- barplot(height = 100*idown["intermediate",,"nonHIV","0.5"], beside = TRUE, 
                    ylab="% reduction in year 10 TB mortality (median [IQR])", 
                    xlab="HIV-unrelated exclusions\nfrom novel regimen", names.arg=c("11%","5%","0%"), cex.lab=1, main="",
                    ylim=c(-25,2), col=cols)
mtext("Novel DR regimen; intermediate efficacy; 5% HIV coprevalence", outer=TRUE, side=3, cex=1.4)
arrows(inoH, 100*idown["intermediate",,"nonHIV","0.25"], inoH, 100*idown["intermediate",,"nonHIV","0.75"], angle=90, code=3, length=0.05)


iH <- barplot(height = 100*idown["intermediate",,"HIV","0.5"], beside = TRUE, 
                  ylab="% reduction in year 10 TB mortality (median [IQR])", 
                  xlab="HIV-related exclusions\nfrom novel regimen", names.arg=c("100%","5%","0%"), cex.lab=1, main="",
                  ylim=c(-25,2), col=cols)
arrows(inoH, 100*idown["intermediate",,"HIV","0.25"], inoH, 100*idown["intermediate",,"HIV","0.75"], angle=90, code=3, length=0.05)


# South Africa exclusions

sexc1 <- read.csv("Exclusions_DRDSTall_SouthAfrica_20160111.idr1.csv")
sexc2 <- read.csv("Exclusions_DRDSTall_SouthAfrica_20160111.idr2.csv")
colnames(sexc2) <- colnames(sexc1)
sexc <- rbind(sexc1, sexc2)
# sdrout <- read.csv("DRcalibration_SouthAfrica_20160111.1.csv"); 
#   sdrout <- sdrout[sdrout[,"rrinc"]/sdrout[,"inc"] > 1/tolerance*sdrout[,"targetdr"] & sdrout[,"rrinc"]/sdrout[,"inc"] < tolerance*sdrout[,"targetdr"], ];
#   drout <- rbind(sdrout[sdrout$idr == 1,], sdrout[sdrout$idr == 2,])
outcome <- "rrdeaths"
sdown <- array(0,dim=c( 3, 3, 2, 5 ));
dimnames(sdown) <- list("efflevel"=levels, "exclevel"=levels, "HIV"=c("HIV","nonHIV"), "q"=c(0.025,0.25,0.5,0.75,0.975))

for (efficacy in levels) for (exclusions in levels) for (H in c("HIV","nonHIV"))
{  sdown[efficacy,exclusions, H,] <- quantile((sexc[ , paste0(outcome, "10", H,"exclusions",exclusions,"efficacy",efficacy)] - 
                                                 sexc[, paste0(outcome,"10allminimal")] )/
                                                     sexc[ , paste0(outcome,"10allminimal")], c(0.025,0.25,0.5,0.75,0.975))
}

par(mar=c(5,5,1,1), mfrow=c(1,2), oma=c(0,0,2,0))

snoH <- barplot(height = 100*sdown["intermediate",,"nonHIV","0.5"], beside = TRUE, 
                ylab="% reduction in year 10 TB mortality (median [IQR])", 
                xlab="HIV-unrelated exclusions\nfrom novel regimen", names.arg=c("11%","5%","0%"), cex.lab=1, main="",
                ylim=c(-40,2), col=cols)
mtext("Novel DR regimen; intermediate efficacy; 60% HIV coprevalence", outer=TRUE, side=3, cex=1.4)
arrows(snoH, 100*sdown["intermediate",,"nonHIV","0.25"], snoH, 100*sdown["intermediate",,"nonHIV","0.75"], angle=90, code=3, length=0.05)


sH <- barplot(height = 100*sdown["intermediate",,"HIV","0.5"], beside = TRUE, 
              ylab="% reduction in year 10 TB mortality (median [IQR])", 
              xlab="HIV-related exclusions\nfrom novel regimen", names.arg=c("100%","5%","0%"), cex.lab=1, main="",
              ylim=c(-40,2), col=cols)
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


#pessimistic adherence baseline scenario

drDST <- numeric(0)
i <- 1; while(file.exists(paste0(location,"DRcalibration_rDSTall.India_20160214p.",i,".csv")))
{drDST <- rbind(drDST, read.csv(paste0(location,"DRcalibration_rDSTall.India_20160214p.",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0
drDST <- drDST[drDST[,"rrinc"]/drDST[,"inc"] > 1/tolerance*drDST[,"targetdr"] & drDST[,"rrinc"]/drDST[,"inc"] < tolerance*drDST[,"targetdr"], ]  #within 3fold if rr incident fraction target
drDST <- drDST[ !is.na(drDST[,"ids"]), ]

novelwide <- numeric(0); i <- 1
while (file.exists(paste0(location,"TRPwideoutput_DSDSTall_rDSTall.India_20160214p.",i,".csv")))
{  novelwide <- rbind(novelwide, read.csv(paste0(location,"TRPwideoutput_DSDSTall_rDSTall.India_20160214p.",i,".csv")))
  i <- i+1 }

outcome <- c("tbdeaths") 

pesdown <- array(0,dim=c( length(elementnames) , 3 , 5 )); 
dimnames(pesdown) <- list("vary"=elementnames, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
for (vary in elementnames) 
{ pesdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drDST[ , paste0(outcome,"10")] )/
                                       drDST[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  pesdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drDST[ , paste0(outcome,"10")] )/
                                       drDST[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
  pesdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drDST[ , paste0(outcome,"10")] )/
                                       drDST[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
}  

par(mar=c(2,4,8,1), mfrow=c(1,1), oma=c(0,0,0,0))
bpctdown <- barplot(height = 100*aperm(pesdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="% reduction (median [IQR])",
                    xlab="", las=2, cex.lab=1, 
                   ylim=c(145*pesdown[1,3,"0.5"], 2), yaxt='n',
                    legend=c("minimal","intermediate","optimal"),
                    args.legend=list(title="Level of varied element(s)", x="bottom",cex=0.9),
                    space=c(0,0,0,1.5,0,0,rep(c(0.5,0,0),5)),
                    col=cols, names.arg=rep("", length(elementnames)))
axis(2, at=seq(-25,0,by=5), labels=paste0(seq(-25,0,by=5),"%"),las=2,cex.axis=0.8)
mtext("Reduction in year 10 TB mortality with novel DS-TB regimens,\n with pessimistic baseline adherence assumptions\nand optimistic 25-50% improvement using novel regimen", cex=1.2, line=4, side=3)
mtext("Varied TRP element(s)", side=3, line=2.5)
text(colMeans(bpctdown) ,0.5, elementlabels, cex=0.8, pos=3, srt=0, font=1, xpd=NA)
arrows(bpctdown, aperm(100*pesdown, c(2,1,3))[,,"0.25"], bpctdown, aperm(100*pesdown, c(2,1,3))[,,"0.75"], angle=90, code=3, length=0.05, xpd=NA)

