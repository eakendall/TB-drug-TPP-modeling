levels <- c("minimal","intermediate","optimal"); 
elementnames <- c("all", set.novelvalues()$elementnames[c(1,4,3,5,2,6)])

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
values <- set.values(); genericvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))

shortelementlabels <- c("All", "% Durably Cured", "Duration", 
                   "Existing Resistance", "Barrier to Resistance", 
                   "Exclusions", "Adherence")[c(1,2,5,4,6,3,7)]

elementlabels <- c("All elements\nvaried", "% Durably Cured\n(optimal conditions)", "Regimen\nDuration", 
                   "Prevalence of\nExisting Resistance\nto Regimen", "Barrier to\nAcquired Novel\nDrug Resistance", 
                   "Medical exclusions,\ncontraindications,\nand early\ndiscontinuations", "Adherence/Burden\nto Patient")[c(1,2,5,4,6,3,7)]
dslabels <- c("","","", "94% cured","97% cured", "99% cured", "6 months","4 months", "2 months", 
              "10% resistant", "3% resistant", "None resistant", "1 acquired resistance per 20 treatments", "1 acquired resistance per 125 treatments", "Minimal acquired resistance",
              "Excludes 100% HIV; 11% non-HIV", "Excludes 10% HIV; 5% non-HIV", "No exclusions", "Same as standard of care", "1.5% fewer dropouts", "3% fewer dropouts")[rep(c(1,4,3,5,2,6), each=3)+c(0,1,2)]
drlabels <- c("","","", "76% cured","94% cured", "97% cured", "20 months","9 months", "6 months", 
              "15% resistant", "5% resistant", "None resistant", "1 acquired resistance per 10 treatments", "1 acquired resistance per 20 treatments", "1 acquired resistance per 125 treatments",
              "Excludes 100% HIV; 11% non-HIV", "Excludes 10% HIV; 5% non-HIV", "No exclusions", "Same as standard of care", "3% fewer dropouts", "6% fewer dropouts")[rep(c(1,4,3,5,2,6), each=3)+c(0,1,2)]
outcomenames <- list("tbdeaths"="% reduction, TB mortality", "rrdeaths" = "% reduction, rifampin-resistant TB mortality", 
                     "panronsets"="Incidence of TB resistant to both all-oral regimens", "nDSTs"= "DSTs performed for novel regimen")
scenarionames <- list("DSDSTall"="with","DSDSTnone"="without", "DRDSTall"="with","DRDSTnone"="without")

longparamnames <- c("Increase in acquired novel drug resistance if resistant to companion drug", 
                    "Probability of also acquiring companion drug resistance if acquiring novel drug resistance",
                    "Increase in poor outcomes if companion-drug resistant", "Increase in poor outcomes if novel-drug resistant", 
                    "Fraction of retreatment patients not receiving rifampin DST", "Fraction of new patients receiving rifampin DST",
                    "Transmission fitness cost of rifampin resistance",
                    "Probability of poor outcome with standard DR regimen",
                    "Probability of acquiring rifampin resistance on standard DS regimen",
                    "Probability of cure by standard DS regimen if rifampin resistant",
                    "Transmission fitness cost of novel-drug resistance",
                    "Increase in relapse probability if completed only 1/3 of treatment course",
                    "Increase in relapse probability if completed only 2/3 of treatment course",
                    "Initial loss to follow up after TB diagnosis", "Rate of loss to follow up during standard regimens",
                    "Rate of active TB diagnosis, new HIV- patients", "Rate of active TB diagnosis, retreatment HIV- patients",
                    "Rate of active TB diagnosis, new HIV+ patients", "Rate of active TB diagnosis, retreatment HIV+ patients",
                    "Fraction of poor outcomes that are relapses (vs failures or deaths)",
                    "Fraction of DS patients not durably cured by standard DS regimen",
                    "Added mortality rate from active TB if HIV-", "Added mortality rate from active TB if HIV+",
                    "Fraction of new TB infections progressing rapidly to active, HIV-", "Fraction of new TB infections progressing rapidly to active, HIV+",
                    "Rate of reactivation from latent TB infection, HIV-", "Rate of reactivation from latent TB infection, HIV+",
                    "Baseline mortality rate", "Added (non-TB) mortality if HIV+",
                    "Reduction in rapid progression of superinfections if prior latent infection",
                    "Rate of relapse to active disease in those relapsing after treatment",
                    "Rate of spontaneous resolution of active TB (HIV- only)",
                    "Transmission coefficient (calibrated to target prevalence)",
                    "HIV infection rate (calibrated to target coprevalence)")
                    


cols <- c("pink", "beige","palegreen")
grays <- c("gray30","gray60","gray90")
blues <- c("blue","lightblue","darkblue")

####### pull in data ###########

alldrout <- numeric(0)
i <- 1; while(file.exists(paste0(location,"DRcalibration_",currenttag,".",i,".csv")))
  {alldrout <- rbind(alldrout, read.csv(paste0(location,"DRcalibration_",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0
alldrout <- alldrout[alldrout[,"rrinc"]/alldrout[,"inc"] > 1/tolerance*alldrout[,"targetdr"] & alldrout[,"rrinc"]/alldrout[,"inc"] < tolerance*alldrout[,"targetdr"], ]  #within 3fold if rr incident fraction target

alldrDST <- numeric(0)
i <- 1; while(file.exists(paste0(location,"DRcalibration_rDSTall.",currenttag,".",i,".csv")))
{alldrDST <- rbind(alldrDST, read.csv(paste0(location,"DRcalibration_rDSTall.",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0
alldrDST <- alldrDST[alldrDST[,"rrinc"]/alldrDST[,"inc"] > 1/tolerance*alldrDST[,"targetdr"] & alldrDST[,"rrinc"]/alldrDST[,"inc"] < tolerance*alldrDST[,"targetdr"], ]  #within 3fold if rr incident fraction target


# trajdrout <- numeric(0)
# i <- 1; while(file.exists(paste0(location,"DRtraj_",currenttag,".",i,".csv")))
# {trajdrout <- rbind(trajdrout, read.csv(paste0(location,"DRtraj_",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0
# trajdrout <- trajdrout[trajdrout[,"X10rrinc"]/trajdrout[,"X10inc"] > 1/tolerance*trajdrout[,"targetdr"] & trajdrout[,"X10rrinc"]/trajdrout[,"X10inc"] < tolerance*trajdrout[,"targetdr"], ]  #within 3fold if rr incident fraction target


allnovelwide <- list()
for (flag in c("DSDSTall_rDSTall.","DSDSTnone_rDSTall.", "DRDSTall_", "DRDSTnone_")) #for (targetepi in names(targetepis))
{
  i <- 1; allnovelwide[[flag]] <- numeric(0)
  while (file.exists(paste0(location,"TRPwideoutput_", flag, currenttag,".",i,".csv")))
  {
    allnovelwide[[flag]] <- rbind(allnovelwide[[flag]], read.csv(paste0(location,"TRPwideoutput_", flag, currenttag,".",i,".csv")))
    i <- i+1
  }
}  


# plotresult <- function(outcome, scenario, elements, barlabels=TRUE, cum=FALSE)
# {
#   
#   if (outcome  ==  "rxtime")
#   { plotrxstacked(outcome, scenario, elements, barlabels, cum)
#   } else if (outcome %in% c("panronsets","dxs"))
#   { plotup(outcome, scenario, elements, barlabels, cum) 
#   } else 
#     { plotpctdown(outcome, scenario, elements, barlabels, cum)
#   }
# }
# 
# # individual plot types:
# 
# plotrxstacked <- function(outcome, scenario, elements, barlabels=TRUE, cum=FALSE)
# {
#   if(missing(novelwide)) novelwide <- allnovelwide[[scenario]]
# 
#   resource <- array(0,dim=c( length(elements) , 3 , 3 )); 
#   dimnames(resource) <- list("vary"=elements, "level"=c("minimal", "intermediate", "optimal"), "reg"=c("First-line","Second-line","Novel"))
#   
#   for (vary in elements) for (nreg in 1:3) 
#   { 
#     outcome <- c("rxtime_s","rxtime_r","rxtime_n", "dxs")[nreg]
#     resource[vary,1,nreg] <- median(novelwide[ , paste0(outcome, "10", vary,"minimal")])
#     resource[vary,2,nreg] <- median(novelwide[ , paste0(outcome, "10allintermediate")])
#     resource[vary,3,nreg] <- median(novelwide[ , paste0(outcome, "10", vary,"optimal")])
#     if (cum==TRUE) {for (t in 1:9)
#     {
#       resource[vary,1,nreg] <- resource[vary,1,nreg] + median(novelwide[ , paste0(outcome, t, vary,"minimal")])
#       resource[vary,2,nreg] <- resource[vary,2,nreg] + median(novelwide[ , paste0(outcome, t,"allintermediate")])
#       resource[vary,3,nreg] <- resource[vary,3,nreg] + median(novelwide[ , paste0(outcome, t, vary,"optimal")])
#     } }
#   }  
#   if (cum==TRUE) {ylab="Cumulative patient-months of treatment\nover 10 years, by regimen"
#   } else  {ylab="Patient-months of treatment in year 10, by regimen"}
#   
#   par(mar=c(7,5,5,1))
#   bres <- barplot(12*aperm(resource,c(1,3,2))[elements,1:3,], beside = FALSE, 
#                   space=c(0.75,0.25,0.25), 
#                   col=grays, ylab=ylab)
#   legend(x = bres[1], y=max(12*rowSums(resource[1,,1:3]))+10, xjust=0.2, yjust=0, fill=rev(grays),
#          c("Novel regimen","Second-line regimen","First-line regimen"), xpd=NA)
#   mtext(paste0("Varying ",elementlabels[which(elements==elementnames)]), side=1, line=5, cex=0.8)
#   
#   return(bres)
# }
# 
# plotup <- function(outcome, scenario, elements, barlabels=TRUE, cum=FALSE, novelwide)
# {
#   if(missing(novelwide)) novelwide <- allnovelwide[[scenario]]
#   
#   up <- array(0,dim=c( length(elements) , 3 , 5)); 
#   dimnames(up) <- list("vary"=elements, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
#   
#   for (vary in elements)
#   { 
#     up[vary,1,] <- quantile( novelwide[ , paste0(outcome, "10", vary,"minimal")], c(0.025,0.25,0.5,0.75,0.975))
#     up[vary,2,] <- quantile( novelwide[ , paste0(outcome, "10allintermediate")] , c(0.025,0.25,0.5,0.75,0.975))
#     up[vary,3,] <- quantile( novelwide[ , paste0(outcome, "10", vary,"optimal")] , c(0.025,0.25,0.5,0.75,0.975))
#     if (cum==TRUE)
#     { for (t in 1:9)
#       {
#         up[vary,1,] <- up[vary,1,] + quantile( novelwide[ , paste0(outcome, "10", vary,"minimal")], c(0.025,0.25,0.5,0.75,0.975))
#         up[vary,2,] <- up[vary,2,] + quantile( novelwide[ , paste0(outcome, "10", vary,"minimal")], c(0.025,0.25,0.5,0.75,0.975))
#         up[vary,3,] <- up[vary,3,] + quantile( novelwide[ , paste0(outcome, "10", vary,"minimal")], c(0.025,0.25,0.5,0.75,0.975))
#       }  
#     }
#   }  
#   
#   par(mar=c(7,4,4,1))
#   bup <- barplot(aperm(up,c(2,1,3))[,,"0.5"], beside = TRUE, 
#                   col=cols, ylab=outcomenames[[outcome]], xlab="", ylim=c(0,max(up[,,"0.75"])), cex.lab=1.2,
#                  main=paste0("Novel ",substr(scenario,1,2)," TB regimen,\n",scenarionames[[scenario]]," DST"))
#   mtext(paste0("Varying ",elementlabels[which(elements==elementnames)]), side=1, line=5, cex=0.8)  
#   
#   arrows(bup, aperm(up, c(2,1,3))[,,"0.25"], bup, aperm(up, c(2,1,3))[,,"0.75"], angle=90, code=3, length=0.05)
#   
#   if (barlabels) 
#   {
#     if (substr(scenario,1,2)=="DS") {text(bup-0.2, max(up[,,"0.75"])/100, dslabels[rep(which(elementnames %in% elements) * 3,each=3) - (2:0)]  ,cex=1, pos=4, srt=90, col="black", font=2)
#     } else if (substr(scenario,1,2)=="DR") {text(bup-0.2, max(up[,,"0.75"])/100, drlabels[rep(which(elementnames %in% elements) * 3,each=3) - (2:0)]  ,cex=1, pos=4, srt=90, col="black", font=2)}
#   }
#   
#   return(bup)
# }  
# 
# plotpctdown <- function(outcome, scenario, elements, barlabels=TRUE, elementlabs=FALSE, main, cum=FALSE, novelwide, drout, drtraj, mar=c(7,4,4,1))
# {
#   par(mar=mar)
#   
#   if (missing(novelwide)) novelwide <- allnovelwide[[scenario]]
#   if (missing(drout)) drout <- alldrout[1:nrow(novelwide),]
#   if (cum==TRUE & missing(drtraj)) drtraj <- trajdrout[1:nrow(novelwide),]
#   
#   if(missing(main)) {main <- paste0("Novel ",substr(scenario,1,2)," TB regimen,\n",scenarionames[[scenario]]," DST")}
#     
#   final_pctdown <- array(0,dim=c( length(elements) , 3 , 5 ));
#   dimnames(final_pctdown) <- list("vary"=elements, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.75,0.975))
#   
#   for (vary in elements) 
#   { 
#     final_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
#                                          drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
#     final_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
#                                          drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
#     final_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
#                                          drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.75,0.975))
#     if(cum==TRUE)
#     {
#       for (t in 1:9)
#       {
#         final_pctdown[vary,1,] <- final_pctdown[vary,1,] + quantile((novelwide[ , paste0(outcome, t, vary,"minimal")] - drout[ , paste0("X",t,".",outcome)] )/
#                                                                       drout[ , paste0("X",t,".",outcome)], c(0.025,0.25,0.5,0.75,0.975))
#         final_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0("X",t,".",outcome)])/
#                                              drout[ , paste0("X",t,".",outcome)], c(0.025,0.25,0.5,0.75,0.975))
#         final_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0("X",t,".",outcome)] )/
#                                              drout[ , paste0("X",t,".",outcome)], c(0.025,0.25,0.5,0.75,0.975))
#         
#         
#       }
#       
#     }
#   }  
# 
#   bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
#                     main=main,
#                     ylim=100*min(final_pctdown[,,"0.5"])* c(1.8,-0.2), ylab=outcomenames[[outcome]], cex.lab=1.2,
#                     col=cols, xlab="" ,names.arg=rep("", length(elements)))
#   arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.25"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.75"], angle=90, code=3, length=0.05)
# 
#   
#   if (barlabels) 
#   {
#     if (substr(scenario,1,2)=="DS") {text(bpctdown+0.3, -0.1, dslabels[rep(which(elementnames %in% elements) * 3,each=3) - (2:0)] ,cex=0.8, pos=2, srt=90, col="black", font=2)
#                                      } else if (substr(scenario,1,2)=="DR") { text(bpctdown+0.4, -0.1, drlabels[rep(which(elementnames %in% elements) * 3,each=3) - (2:0)]  ,cex=1, pos=2, srt=90, col="black", font=2)  }
#   }
#   if (elementlabs)
#   {  text(colMeans(bpctdown)-0.5 ,-min(final_pctdown[,,"0.5"]), elementlabels, cex=0.8, pos=4, srt=90, font=1, xpd=NA) }
#   
#   
#   if(length(elements)==1) mtext(paste0("Varying ",elementlabels[which(elements==elementnames)]), side=1, cex=0.8, line=5)
# 
#   
#   return(list("positions"=bpctdown,"array"=final_pctdown))
# }  