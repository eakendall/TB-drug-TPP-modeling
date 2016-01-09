source("TPPmat.R")
currenttag <- "India_20160105"; tolerance <- 1.5; location=""
levels <- c("minimal","intermediate","optimal"); 
elementnames <- c("all", set.novelvalues()$elementnames)

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]

elementlabels <- c("All elements\nvaried", "% Durably Cured\n(optimal conditions)", "Regimen Duration", 
                   "Prevalence of\nExisting Resistance\nto Regimen", "Barrier to\nAcquired Novel\nDrug Resistance", 
                   "Exclusions and\nContraindications", "Adherence/Burden\nto Patient")
dslabels <- c("","","", "94% cured","97% cured", "99% cured", "6 months","4 months", "2 months", 
              "10% resistant", "3% resistant", "None resistant", "1 acquired resistance per 20 treatments", "1 acquired resistance per 125 treatments", "Minimal acquired resistance",
              "Excludes 100% HIV; 11% non-HIV", "Excludes 10% HIV; 5% non-HIV", "No exclusions", "Same burden as standard of care", "1.5% fewer dropouts", "3% fewer dropouts")
drlabels <- c("","","", "76% cured","94% cured", "97% cured", "20 months","9 months", "6 months", 
              "15% resistant", "5% resistant", "None resistant", "1 acquired resistance per 10 treatments", "1 acquired resistance per 20 treatments", "1 acquired resistance per 125 treatments",
              "Excludes 100% HIV; 11% non-HIV", "Excludes 10% HIV; 5% non-HIV", "No exclusions", "Same burden as standard of care", "3% fewer dropouts", "6% fewer dropouts")
outcomenames <- list("tbdeaths"="% reduction, TB mortality", "rrdeaths" = "% reductions, rifampin-resistant TB mortality", 
                     "panronsets"="Incidence of TB untreatable with all-oral regimens", "nDSTs"= "DSTs performed for novel regimen")
scenarionames <- list("DSDSTall"="with","DSDSTnone"="without", "DRDSTall"="with","DRSDSTnone"="without")
cols <- c("pink", "beige","palegreen")
grays <- c("gray30","gray60","gray90")

####### pull in data ###########

alldrout <- numeric(0)
i <- 1; while(file.exists(paste0(location,"DRcalibration_",currenttag,".",i,".csv")))
  {alldrout <- rbind(alldrout, read.csv(paste0(location,"DRcalibration_",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0
alldrout <- alldrout[alldrout[,"rrinc"]/alldrout[,"inc"] > 1/tolerance*alldrout[,"targetdr"] & alldrout[,"rrinc"]/alldrout[,"inc"] < tolerance*alldrout[,"targetdr"], ]  #within 3fold if rr incident fraction target


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

# outcomes <- rep(c("tbdeaths", "panronsets", "rrdeaths", "panronsets"), each=2); 
outcomes <- rep(c("tbdeaths","rrdeaths", "panronsets", "rxtime"), 2)
# scenarios <- c("DSDSTall", "DSDSTnone","DSDSTall", "DSDSTnone","DRDSTall","DRDSTnone","DRDSTall","DRDSTnone"); 
scenarios <- c("DSDSTall", "DSDSTall","DSDSTall", "DSDSTall","DRDSTall","DRDSTall","DRDSTall","DRDSTall"); #wil use the above instead
plotelements <- rep("barrier",8)
par(mfrow=c(2,4))
for (version in 1:8) plotresult(outcomes[version], scenarios[version], plotelements[version])

plotresult <- function(outcome, scenario, elements, barlabels=TRUE, cum=FALSE)
{
  
  if (outcome  ==  "rxtime")
  { plotrxstacked(outcome, scenario, elements, barlabels, cum)
  } else if (outcome %in% c("panronsets","dxs"))
  { plotup(outcome, scenario, elements, barlabels, cum) 
  } else 
    { plotpctdown(outcome, scenario, elements, barlabels, cum)
  }
}

# individual plot types:

plotrxstacked <- function(outcome, scenario, elements, barlabels=TRUE, cum=FALSE)
{
  novelwide <- trpwide[[scenario]]

  resource <- array(0,dim=c( length(elements) , 3 , 3 )); 
  dimnames(resource) <- list("vary"=elements, "level"=c("minimal", "intermediate", "optimal"), "reg"=c("First-line","Second-line","Novel"))
  
  for (vary in elements) for (nreg in 1:3) 
  { 
    outcome <- c("rxtime_s","rxtime_r","rxtime_n", "dxs")[nreg]
    resource[vary,1,nreg] <- median(novelwide[ , paste0(outcome, "10", vary,"minimal")])
    resource[vary,2,nreg] <- median(novelwide[ , paste0(outcome, "10allintermediate")])
    resource[vary,3,nreg] <- median(novelwide[ , paste0(outcome, "10", vary,"optimal")])
    if (cum==TRUE) {for (t in 1:9)
    {
      resource[vary,1,nreg] <- resource[vary,1,nreg] + median(novelwide[ , paste0(outcome, t, vary,"minimal")])
      resource[vary,2,nreg] <- resource[vary,2,nreg] + median(novelwide[ , paste0(outcome, t,"allintermediate")])
      resource[vary,3,nreg] <- resource[vary,3,nreg] + median(novelwide[ , paste0(outcome, t, vary,"optimal")])
    } }
  }  
  if (cum==TRUE) {ylab="Cumulative patient-months of treatment\nover 10 years, by regimen"
  } else  {ylab="Patient-months of treatment in year 10, by regimen"}
  
  par(mar=c(7,5,5,1))
  bres <- barplot(12*aperm(resource,c(1,3,2))[elements,1:3,], beside = FALSE, 
                  space=c(0.75,0.25,0.25), 
                  col=grays, ylab=ylab)
  legend(x = bres[1], y=max(12*rowSums(resource[1,,1:3]))+10, xjust=0.2, yjust=0, fill=rev(grays),
         c("Novel regimen","Second-line regimen","First-line regimen"), xpd=NA)
  mtext(paste0("Varying ",elementlabels[which(elements==elementnames)]), side=1, line=5, cex=0.8)
}

plotup <- function(outcome, scenario, elements, barlabels=TRUE, cum=FALSE)
{
  novelwide <- trpwide[[scenario]]
  
  up <- array(0,dim=c( length(elements) , 3 , 5)); 
  dimnames(up) <- list("vary"=elements, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))
  
  for (vary in elements)
  { 
    up[vary,1,] <- quantile( novelwide[ , paste0(outcome, "10", vary,"minimal")], c(0.025,0.25,0.5,0.075,0.975))
    up[vary,2,] <- quantile( novelwide[ , paste0(outcome, "10allintermediate")] , c(0.025,0.25,0.5,0.075,0.975))
    up[vary,3,] <- quantile( novelwide[ , paste0(outcome, "10", vary,"optimal")] , c(0.025,0.25,0.5,0.075,0.975))
    if (cum==TRUE)
    { for (t in 1:9)
      {
        up[vary,1,] <- up[vary,1,] + quantile( novelwide[ , paste0(outcome, "10", vary,"minimal")], c(0.025,0.25,0.5,0.075,0.975))
        up[vary,2,] <- up[vary,2,] + quantile( novelwide[ , paste0(outcome, "10", vary,"minimal")], c(0.025,0.25,0.5,0.075,0.975))
        up[vary,3,] <- up[vary,3,] + quantile( novelwide[ , paste0(outcome, "10", vary,"minimal")], c(0.025,0.25,0.5,0.075,0.975))
      }  
    }
  }  

  bup <- barplot(aperm(up,c(2,1,3))[,,"0.5"], beside = TRUE, 
                  col=cols, ylab=outcomenames[[outcome]], xlab="", ylim=c(0,max(up[,,"0.975"])), cex.lab=1.4,
                 main=paste0("Novel ",substr(scenario,1,2)," TB regimen,\n",scenarionames[[scenario]]," DST"))
  mtext(paste0("Varying ",elementlabels[which(elements==elementnames)]), side=1, line=5, cex=0.8)  
  
  arrows(bup, aperm(up, c(2,1,3))[,,"0.025"], bup, aperm(up, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)
  
  if (barlabels) 
  {
    if (substr(scenario,1,2)=="DS") {text(bup-0.2, max(up[,,"0.975"])/100, dslabels[rep(which(elementnames %in% elements) * 3,each=3) - (2:0)]  ,cex=1, pos=4, srt=90, col="black", font=2)
    } else if (substr(scenario,1,2)=="DR") {text(bup-0.2, max(up[,,"0.975"])/100, drlabels[rep(which(elementnames %in% elements) * 3,each=3) - (2:0)]  ,cex=1, pos=4, srt=90, col="black", font=2)}
  }
}  

plotpctdown <- function(outcome, scenario, elements, barlabels=TRUE, cum=FALSE, novelwide, drout)
{
  if (missing(novelwide)) novelwide <- trpwide[[scenario]]
  if (missing(drout)) drout <- alldrout[1:nrow(novelwide),]
  
  final_pctdown <- array(0,dim=c( length(elements) , 3 , 5 ));
  dimnames(final_pctdown) <- list("vary"=elements, "level"=c("minimal", "intermediate", "optimal"), "q"=c(0.025,0.25,0.5,0.075,0.975))
  
  for (vary in elements) 
  { 
    final_pctdown[vary,1,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"minimal")] - drout[ , paste0(outcome,"10")] )/
                                         drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
    final_pctdown[vary,2,] <- quantile((novelwide[ , paste0(outcome, "10allintermediate")] - drout[ , paste0(outcome,"10")] )/
                                         drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
    final_pctdown[vary,3,] <- quantile((novelwide[ , paste0(outcome, "10", vary,"optimal")] - drout[ , paste0(outcome,"10")] )/
                                         drout[ , paste0(outcome,"10")], c(0.025,0.25,0.5,0.075,0.975))
  }  

  bpctdown <- barplot(height = 100*aperm(final_pctdown, c(2,1,3))[,,"0.5"], beside = TRUE, 
                    main=paste0("Novel ",substr(scenario,1,2)," TB regimen,\n",scenarionames[[scenario]]," DST"),
                    ylim=100*min(final_pctdown[,,"0.5"])* c(1.5,-0.2), ylab=outcomenames[[outcome]],
                    col=cols, xlab="")
  arrows(bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.025"], bpctdown, aperm(100*final_pctdown, c(2,1,3))[,,"0.975"], angle=90, code=3, length=0.05)

  
  if (barlabels) 
  {
    if (substr(scenario,1,2)=="DS") {text(bpctdown+0.4, -0.1, dslabels[rep(which(elementnames %in% elements) * 3,each=3) - (2:0)] ,cex=1, pos=2, srt=90, col="black", font=2)
                                     } else if (substr(scenario,1,2)=="DR") { text(bpctdown+0.4, -0.1, drlabels[rep(which(elementnames %in% elements) * 3,each=3) - (2:0)]  ,cex=1, pos=2, srt=90, col="black", font=2)  }
  }
  if(length(elements)==1) mtext(paste0("Varying ",elementlabels[which(elements==elementnames)]), side=1, cex=0.8, line=5)
}  