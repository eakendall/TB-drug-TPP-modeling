currenttag <- "20151222b"

levels <- c("minimal","intermediate","optimal"); 
outcome <- c("tbdeaths") #can set up loop over multiple outcomes
elementnames <- c("all", set.novelvalues()$elementnames)

drout <- screendrout()

novelwide <- list()
for (targetpt in c("DS","DR")) for (DST in c("DSTall","DSTnone")) for (targetepi in names(targetepis))
{
  novewide[[paste0(targetpt, DST)]] <- loadnovel(targetpt=targetpt, DST=DST, targetepi=targetepi, tag=currenttag)$wide
  
  

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

plot(0:10, traj[,"minimal","0.5"], ylim=c(0,100))
points(0:10, traj[,"intermediate","0.5"])
points(0:10, traj[,"optimal","0.5"])
points(10,median(drout[,"tbdeaths10"]))

par(mar=c(7,4,3,3))
bdiff <- barplot(height = aperm(final_diff, c(2,1,3))[,,"0.5"], beside = TRUE, ylab="change in TB mortality compared to typical TRP", xlab="", las=2, ylim=c(-6,6),
        legend=c("minimal","optimal"), args.legend=list(title="Level of varied element(s)",x="topright",cex=0.8)); mtext("Varied TRP element(s)", side=1, line=6)
arrows(bdiff, aperm(final_diff, c(2,1,3))[,,"0.025"], bdiff, aperm(final_diff, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)


par(mar=c(7,5,3,3))
bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    ylab="Percent reduction in TB mortality with novel regimen\n compared to projection under continued current practice", xlab="", las=2, cex.lab=0.9, 
                 legend=c("minimal","intermediate","optimal"), args.legend=list(title="Level of varied element(s)",x="bottomright",cex=0.8),
                 main=paste0("Novel regimen for DS TB, with universal DST"), cex.main=0.9); 
mtext("Varied TRP element(s)", side=1, line=6)
arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)

# Make this 2x2 for DS/DR and w/wo DST


# resource use: barplots as above but for outcomes of diagnoses, DSTs (rif and novel in same plot), and rxmonths (all 3 in same plot)
  
