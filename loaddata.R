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
                   "Pre-Existing Resistance", "Barrier to Resistance", 
                   "Medical Contraindications", "Tolerability/Adherence","Scalability","Increase in RR diagnosis")[c(1,2,5,4,6,3,7,8,9)]

elementlabels <- c("All elements\nvaried\nincluding scaleup", "Efficacy\n(% durably cured,\nif susceptible and\ncomplete treatment)", "Regimen\nDuration", 
                   "Prevalence of\nExisting Resistance\nto Regimen", "Barrier to\nAcquired Novel\nDrug Resistance", 
                   "Medical exclusions,\ncontraindications,\nand adverse reactions", "Adherence/\nTolerability", 
                   "Uptake of novel\nregimen scale-up\n(to replace\nSOC regimen)", "Associated\nscale-up of\nRR diagnosis\nand treatment")[c(1,2,5,4,6,3,7,8,9)]
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
                    "Fraction of RR patients not durably cured by standard RR regimen",
                    "Probability of acquiring rifampin resistance on standard RS regimen",
                    "Probability of cure by standard RS regimen if rifampin resistant",
                    "Time until new treatment, for patients failing current regimen",
                    "Transmission fitness cost of novel-drug resistance",
                    "Increase in relapse probability if completed only 1/3 of treatment course",
                    "Increase in relapse probability if completed only 2/3 of treatment course",
                    "Initial loss to follow up after TB diagnosis", "Rate of loss to follow up during standard regimens",
                    "Rate of active TB diagnosis, new HIV- patients", "Rate of active TB diagnosis, retreatment HIV- patients",
                    "Rate of active TB diagnosis, new HIV+ patients", "Rate of active TB diagnosis, retreatment HIV+ patients",
                    "Fraction of poor outcomes that are relapses (vs failures or deaths)",
                    "Fraction of RS patients not durably cured by standard RS regimen",
                    "Added mortality rate from active TB (HIV-)", "Added mortality rate from active TB (HIV+)",
                    "Fraction of new TB infections progressing rapidly to active, HIV-", "Fraction of new TB infections progressing rapidly to active, HIV+",
                    "Rate of reactivation from latent TB infection, HIV-", "Rate of reactivation from latent TB infection, HIV+",
                    "Baseline mortality rate", "Added (non-TB) mortality if HIV+",
                    "Reduction in rapid progression of superinfections if prior latent infection",
                    "Rate of relapse to active disease in those relapsing after treatment",
                    "Rate of spontaneous resolution of active TB (HIV- only)",
                    "Transmission coefficient(calibrated to target prevalence)",
                    "HIV infection rate (calibrated to target coprevalence)")
                    
paramsigns <- c(1,
                1,
                1,1,
                -1,1,
                -1,
                -1,
                1,
                1,
                1,
                -1,
                1,
                1,
                1,1,
                1,1,
                1,1,
                1,
                -1,
                1,1,
                1,1,
                1,1,
                1,1,
                1,
                1,
                1,
                1,
                1)
brokenlongparamnames <- c("Increase in acquired novel drug resistance if resistant to companion drug", 
                    "Probability of also acquiring companion drug resistance if acquiring novel drug resistance",
                    "Increase in poor outcomes if companion-drug resistant", "Increase in poor outcomes if novel-drug resistant", 
                    "Fraction of RR detected,\nretreatment patients", "Fraction of RR detected,\nnew patients",
                    "Transmissibility of\nresistant strains",
                    "Efficacy of standard\nRR regimen",
                    "Probability of acquiring RR\non standard RS regimen",
                    "Efficacy of standard RS\nregimen against RR TB",
                    "Time to new treatment after\nfailing current regimen",
                    "Transmissibility of novel-drug\nresistant strains",
                    "Extra relapse from nonadherence,\n if 1/3 of treatment taken",
                    "Extra relapse from nonadherence,\n if 2/3 of treatment taken",
                    "Pre-treatment loss to follow up", "Monthly loss to follow up,\nstandard regimens",
                    "New TB diagnosis rate, HIV-", "TB re-diagnosis rate, HIV-",
                    "TB diagnosis rate,\nnew HIV+ patients", "TB diagnosis rate,\nretreatment HIV+ patients",
                    "Relapse as proportion of\nthose not cured",
                    "Efficacy of standard\nRS regimen",
                    "TB mortality rate, HIV-", "TB mortality rate, HIV+",
                    "Fraction progressing\nrapidly, HIV-", "Fraction progressing\nrapidly, HIV+",
                    "TB reactivation rate, HIV-", "TB reactivation, HIV+",
                    "Baseline mortality rate, HIV-", "Baseline mortality rate, HIV+",
                    "Reduction in rapid progression\nif had prior latent infection",
                    "Rate of relapse after\nineffective treatment",
                    "Rate of spontaneous\nresolution (HIV- only)",
                    "Transmission coefficient (calibrated\nto target TB prevalence)",
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

drtraj_ds <- numeric(0)
i <- 1; while(file.exists(paste0(location,"DRtraj_rDSTall.",currenttag,".",i,".csv")))
{drtraj_ds <- rbind(drtraj_ds, read.csv(paste0(location,"DRtraj_rDSTall.",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0
drtraj_ds <- drtraj_ds[drtraj_ds[,"X0rrinc"]/drtraj_ds[,"X0inc"] > 1/tolerance*drtraj_ds[,"targetdr"] & drtraj_ds[,"X0rrinc"]/drtraj_ds[,"X0inc"] < tolerance*drtraj_ds[,"targetdr"], ]  #within 3fold if rr incident fraction target

drtraj_dr <- numeric(0)
i <- 1; while(file.exists(paste0(location,"DRtraj_",currenttag,".",i,".csv")))
{drtraj_dr <- rbind(drtraj_dr, read.csv(paste0(location,"DRtraj_",currenttag,".",i,".csv"), header = TRUE)); i <- i+1} #saved results from dr sampling runs at time 0
drtraj_dr <- drtraj_dr[drtraj_dr[,"X0rrinc"]/drtraj_dr[,"X0inc"] > 1/tolerance*drtraj_dr[,"targetdr"] & drtraj_dr[,"X0rrinc"]/drtraj_dr[,"X0inc"] < tolerance*drtraj_dr[,"targetdr"], ]  #within 3fold if rr incident fraction target



allnovelwide <- list()
for (flag in c("DSDSTall_rDSTall.","DRDSTall_"))#, "DSDSTnone_rDSTall.", "DRDSTnone_")) 
{
  i <- 1; allnovelwide[[flag]] <- numeric(0)
  while (file.exists(paste0(location,"TRPwideoutput_", flag, currenttag,".",i,".csv")))
  {
    allnovelwide[[flag]] <- rbind(allnovelwide[[flag]], read.csv(paste0(location,"TRPwideoutput_", flag, currenttag,".",i,".csv")))
    i <- i+1
  }
  allnovelwide[[flag]] <- allnovelwide[[flag]][ !is.na(allnovelwide[[flag]][,"targetprev"]), ]
}  
novelwide1 <- allnovelwide[["DSDSTall_rDSTall."]]
novelwide3 <- allnovelwide[["DRDSTall_"]]
novelwide3 <- novelwide3[order(novelwide3$ids, novelwide3$idr),]
rm(allnovelwide)

only_ds <- numeric(0); i<- 1 
while(file.exists(paste0("Only_DSDSTall_",currenttag,".",i,".csv"))) { only_ds <- rbind(only_ds, read.csv(paste0("Only_DSDSTall_",currenttag,".",i,".csv"))); i <- i+1 }

intonly_ds <- numeric(0); i<- 1 
while(file.exists(paste0("IntOnly_DSDSTall_",currenttag,".",i,".csv"))) { intonly_ds <- rbind(intonly_ds, read.csv(paste0("IntOnly_DSDSTall_",currenttag,".",i,".csv"))); i <- i+1 }

allbut_ds <- numeric(0); i<- 1 
while(file.exists(paste0("Allbut_DSDSTall_",currenttag,".",i,".csv"))) { allbut_ds <- rbind(allbut_ds, read.csv(paste0("Allbut_DSDSTall_",currenttag,".",i,".csv"))); i <- i+1 }

intallbut_ds <- numeric(0); i<- 1 
while(file.exists(paste0("IntAllbut_DSDSTall_",currenttag,".",i,".csv"))) { intallbut_ds <- rbind(intallbut_ds, read.csv(paste0("IntAllbut_DSDSTall_",currenttag,".",i,".csv"))); i <- i+1 }

only_dr <- numeric(0); i<- 1 
while(file.exists(paste0("Only_DRDSTall_",currenttag,".",i,".csv"))) { only_dr <- rbind(only_dr, read.csv(paste0("Only_DRDSTall_",currenttag,".",i,".csv"))); i <- i+1 }

only_dr2 <- only_dr[!duplicated(1000*only_dr$ids + only_dr$idr),] 
only_dr <- only_dr2

intonly_dr <- numeric(0); i<- 1 
while(file.exists(paste0("IntOnly_DRDSTall_",currenttag,".",i,".csv"))) { intonly_dr <- rbind(intonly_dr, read.csv(paste0("IntOnly_DRDSTall_",currenttag,".",i,".csv"))); i <- i+1 }

intonly_dr2 <- intonly_dr[!duplicated(1000*intonly_dr$ids + intonly_dr$idr),] 
intonly_dr <- intonly_dr2

allbut_dr <- numeric(0); i<- 1 
while(file.exists(paste0("Allbut_DRDSTall_",currenttag,".",i,".csv"))) { allbut_dr <- rbind(allbut_dr, read.csv(paste0("Allbut_DRDSTall_",currenttag,".",i,".csv"))); i <- i+1 }
allbut_dr <- allbut_dr[!duplicated(allbut_dr[,c("ids","idr")]),]

intallbut_dr <- numeric(0); i<- 1 
while(file.exists(paste0("IntAllbut_DRDSTall_",currenttag,".",i,".csv"))) { intallbut_dr <- rbind(intallbut_dr, read.csv(paste0("IntAllbut_DRDSTall_",currenttag,".",i,".csv"))); i <- i+1 }


allminopt_ds <- numeric(0); i <- 1
while(file.exists(paste0("Allminopt_DSDSTall_",currenttag,".",i,".csv"))) { allminopt_ds <- rbind(allminopt_ds, read.csv(paste0("Allminopt_DSDSTall_",currenttag,".",i,".csv"))); i <- i+1 }

allminopt_dr <- numeric(0); i <- 1
while(file.exists(paste0("Allminopt_DRDSTall_",currenttag,".",i,".csv"))) { allminopt_dr <- rbind(allminopt_dr, read.csv(paste0("Allminopt_DRDSTall_",currenttag,".",i,".csv"))); i <- i+1 }

# varied all here instead of only scaleup
sens_ds <- numeric(0); i <- 1
while(file.exists(paste0("Scaleupminopt_DSDSTall_",currenttag,".",i,".csv"))) { sens_ds <- rbind(sens_ds, read.csv(paste0("Scaleupminopt_DSDSTall_",currenttag,".",i,".csv"))); i <- i+1 }


allbut_hiv <- numeric(0); i <- 1
while(file.exists(paste0("HIVAllbut_DSDSTall_",currenttag,".",i,".csv"))) { allbut_hiv <- rbind(allbut_hiv, read.csv(paste0("HIVAllbut_DSDSTall_",currenttag,".",i,".csv"))); i <- i+1 }

combos_dr <- numeric(0); i <- 1
while(file.exists(paste0("TRPcombos_DRDSTall_",currenttag,".",i,".csv"))) { combos_dr <- rbind(combos_dr, read.csv(paste0("TRPcombos_DRDSTall_",currenttag,".",i,".csv"))); i <- i+1 }
combos_dr <- combos_dr[!duplicated(combos_dr[,c("ids","idr")]),]

sens_dr <- numeric(0); i <- 1
while(file.exists(paste0("Scaleupminopt_DRDSTall_",currenttag,".",i,".csv"))) { sens_dr <- rbind(sens_dr, read.csv(paste0("Scaleupminopt_DRDSTall_",currenttag,".",i,".csv"))); i <- i+1 }
