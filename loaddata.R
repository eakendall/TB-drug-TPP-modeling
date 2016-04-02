levels <- c("minimal","intermediate","optimal"); 
elementnames <- c("all", set.novelvalues()$elementnames[c(1,4,3,5,2,6,7,8)])
dselementnames <- c("all", set.novelvalues()$elementnames[c(1,4,3,5,2,6,7)])
drelementnames <- c("all", set.novelvalues()$elementnames[c(1,4,3,5,2,6,8)])

dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE)
drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE)
novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
values <- set.values(); genericvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))

shortelementlabels <- c("All", "Efficacy", "Duration", 
                   "Existing Resistance", "Barrier to Resistance", 
                   "Exclusions", "Adherence","Scalability","Increase in RR diagnosis")[c(1,2,5,4,6,3,7,8,9)]

elementlabels <- c("All elements\nvaried", "Efficacy\n(% durably cured,\nif susceptible and\ncomplete treatment)", "Regimen\nDuration", 
                   "Prevalence of\nExisting Resistance\nto Regimen", "Barrier to\nAcquired Novel\nDrug Resistance", 
                   "Medical exclusions,\ncontraindications,\nand adverse reactions", "Adherence/\nTolerability", 
                   "Reach of novel\nregimen scale-up\n(to replace\nSOC regimen)", "Associated\nexpansion of\nRR diagnosis\nand treatment")[c(1,2,5,4,6,3,7,8,9)]
dselementlabels <- elementlabels[-9]
drelementlabels <- elementlabels[-8]
dslabels <- c("","","", "94% cured","97% cured", "99% cured", "6 months","4 months", "2 months", 
              "10% resistant", "3% resistant", "None resistant", "1 acquired resistance per 20 treatments", "1 acquired resistance per 125 treatments", "Minimal acquired resistance",
              "Excludes 100% HIV; 11% non-HIV", "Excludes 10% HIV; 5% non-HIV", "No exclusions", "Same as standard of care", "1.5% fewer dropouts", "3% fewer dropouts",
              "Reaches 50% of eligible treated patients", "Reaches 75% of eligible treated patients", "Reaches 100% of eligible treated patients", 
              "","","")[rep(c(1,4,3,5,2,6,7,8), each=3)+c(0,1,2)]
drlabels <- c("","","", "76% cured","94% cured", "97% cured", "20 months","9 months", "6 months", 
              "15% resistant", "5% resistant", "None resistant", "1 acquired resistance per 10 treatments", "1 acquired resistance per 20 treatments", "1 acquired resistance per 125 treatments",
              "Excludes 100% HIV; 11% non-HIV", "Excludes 10% HIV; 5% non-HIV", "No exclusions", "Same as standard of care", "3% fewer dropouts", "6% fewer dropouts",
              "Reaches 50% of eligible treated patients", "Reaches 75% of eligible treated patients", "Reaches 100% of eligible treated patients",
              "No increase in RR diagnosis", "Allows 50% reduction in RR under-diagnosis", "Allows universal RR diagnosis")[rep(c(1,4,3,5,2,6), each=3)+c(0,1,2)]
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

# alldrout <- numeric(0)
# i <- 1; while(file.exists(paste0(location,"DRcalibration_",currenttag,".",i,".csv")))
#   {alldrout <- rbind(alldrout, read.csv(paste0(location,"DRcalibration_",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0
# alldrout <- alldrout[alldrout[,"rrinc"]/alldrout[,"inc"] > 1/tolerance*alldrout[,"targetdr"] & alldrout[,"rrinc"]/alldrout[,"inc"] < tolerance*alldrout[,"targetdr"], ]  #within 3fold if rr incident fraction target
# 
# alldrDST <- numeric(0)
# i <- 1; while(file.exists(paste0(location,"DRcalibration_rDSTall.",currenttag,".",i,".csv")))
# {alldrDST <- rbind(alldrDST, read.csv(paste0(location,"DRcalibration_rDSTall.",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0
# alldrDST <- alldrDST[alldrDST[,"rrinc"]/alldrDST[,"inc"] > 1/tolerance*alldrDST[,"targetdr"] & alldrDST[,"rrinc"]/alldrDST[,"inc"] < tolerance*alldrDST[,"targetdr"], ]  #within 3fold if rr incident fraction target

redo_alldrDST <- numeric(0)
i <- 1; while(file.exists(paste0(location,"DRcalibration_redo_rDSTall.",currenttag,".",i,".csv")))
{redo_alldrDST <- rbind(redo_alldrDST, read.csv(paste0(location,"DRcalibration_redo_rDSTall.",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0
redo_drDST <- redo_alldrDST[redo_alldrDST[,"rrinc"]/redo_alldrDST[,"inc"] > 1/tolerance*redo_alldrDST[,"targetdr"] & redo_alldrDST[,"rrinc"]/redo_alldrDST[,"inc"] < tolerance*redo_alldrDST[,"targetdr"], ]  #within 3fold if rr incident fraction target

redo_alldrout <- numeric(0)
i <- 1; while(file.exists(paste0(location,"DRcalibration_redo_",currenttag,".",i,".csv")))
{redo_alldrout <- rbind(redo_alldrout, read.csv(paste0(location,"DRcalibration_redo_",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0
redo_drout <- redo_alldrout[redo_alldrout[,"rrinc"]/redo_alldrout[,"inc"] > 1/tolerance*redo_alldrout[,"targetdr"] & redo_alldrout[,"rrinc"]/redo_alldrout[,"inc"] < tolerance*redo_alldrout[,"targetdr"], ]  #within 3fold if rr incident fraction target

alldrout <- redo_drout
alldrDST <- redo_drDST


allnovelwide <- list()
for (flag in c("DSDSTall_rDSTall.","DRDSTall_", "DSDSTnone_rDSTall.", "DRDSTnone_")) #for (targetepi in names(targetepis))
{
  i <- 1; allnovelwide[[flag]] <- numeric(0)
  while (file.exists(paste0(location,"TRPwideoutput_", flag, currenttag,".",i,".csv")))
  {
    allnovelwide[[flag]] <- rbind(allnovelwide[[flag]], read.csv(paste0(location,"TRPwideoutput_", flag, currenttag,".",i,".csv")))
    i <- i+1
  }
  allnovelwide[[flag]] <- allnovelwide[[flag]][ !is.na(allnovelwide[[flag]][,"targetprev"]), ]
}  


only_ds <- numeric(0); i<- 1 
while(file.exists(paste0("Only_DSDSTall_India_20160313p.",i,".csv"))) { only_ds <- rbind(only_ds, read.csv(paste0("Only_DSDSTall_India_20160313p.",i,".csv"))); i <- i+1 }

intonly_ds <- numeric(0); i<- 1 
while(file.exists(paste0("IntOnly_DSDSTall_India_20160313p.",i,".csv"))) { intonly_ds <- rbind(intonly_ds, read.csv(paste0("IntOnly_DSDSTall_India_20160313p.",i,".csv"))); i <- i+1 }

allbut_ds <- numeric(0); i<- 1 
while(file.exists(paste0("Allbut_DSDSTall_India_20160313p.",i,".csv"))) { allbut_ds <- rbind(allbut_ds, read.csv(paste0("Allbut_DSDSTall_India_20160313p.",i,".csv"))); i <- i+1 }

intallbut_ds <- numeric(0); i<- 1 
while(file.exists(paste0("IntAllbut_DSDSTall_India_20160313p.",i,".csv"))) { intallbut_ds <- rbind(intallbut_ds, read.csv(paste0("IntAllbut_DSDSTall_India_20160313p.",i,".csv"))); i <- i+1 }

only_dr <- numeric(0); i<- 1 
while(file.exists(paste0("Only_DRDSTall_India_20160313p.",i,".csv"))) { only_dr <- rbind(only_dr, read.csv(paste0("Only_DRDSTall_India_20160313p.",i,".csv"))); i <- i+1 }

intonly_dr <- numeric(0); i<- 1 
while(file.exists(paste0("IntOnly_DRDSTall_India_20160313p.",i,".csv"))) { intonly_dr <- rbind(intonly_dr, read.csv(paste0("IntOnly_DRDSTall_India_20160313p.",i,".csv"))); i <- i+1 }

allbut_dr <- numeric(0); i<- 1 
while(file.exists(paste0("Allbut_DRDSTall_India_20160313p.",i,".csv"))) { allbut_dr <- rbind(allbut_dr, read.csv(paste0("Allbut_DRDSTall_India_20160313p.",i,".csv"))); i <- i+1 }

intallbut_dr <- numeric(0); i<- 1 
while(file.exists(paste0("IntAllbut_DRDSTall_India_20160313p.",i,".csv"))) { intallbut_dr <- rbind(intallbut_dr, read.csv(paste0("IntAllbut_DRDSTall_India_20160313p.",i,".csv"))); i <- i+1 }


allminopt_ds <- numeric(0); i <- 1
while(file.exists(paste0("Allminopt_DSDSTall_India_20160313p.",i,".csv"))) { allminopt_ds <- rbind(allminopt_ds, read.csv(paste0("Allminopt_DSDSTall_India_20160313p.",i,".csv"))); i <- i+1 }

allminopt_dr <- numeric(0); i <- 1
while(file.exists(paste0("Allminopt_DRDSTall_India_20160313p.",i,".csv"))) { allminopt_dr <- rbind(allminopt_dr, read.csv(paste0("Allminopt_DRDSTall_India_20160313p.",i,".csv"))); i <- i+1 }
