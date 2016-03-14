targetpt <- "DS"; DST <- "DSTall"
drout <- read.csv(paste0(location,"DRcalibration_rDSTall.India_20160309p.1.csv")); drrow <- 1


source("TPPmat.R")
dssetup <- setup.model(DRera=FALSE, treatSL=FALSE, treatnovel=FALSE); drsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=FALSE); novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
values <- set.values(); mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]; elementnames <- c("all", set.novelvalues()$elementnames); if (targetpt=="DS") elementnames <- elementnames[-which(elementnames=="riftest")]

tolerance <- 1.5; drout <- drout[drout[,"rrinc"]/drout[,"inc"] > 1/tolerance*drout[,"targetdr"] & drout[,"rrinc"]/drout[,"inc"] < tolerance*drout[,"targetdr"], ]  #within 3fold if rr incident fraction target
valuevect <- unlist(drout[drrow, 5+(1:length(unlist(mergedvalues)))]) ; names(valuevect) <- names(unlist(mergedvalues))
v <- 0; for (pname in names(mergedvalues))
{ mergedvalues[[pname]] <- valuevect[v+(1:length(mergedvalues[[pname]]))]
  v <- v + length(mergedvalues[[pname]]) } # move forward to start of next par vector in sampled values }    

values <- sampleTRP(mergedvalues = mergedvalues, targetpt=targetpt, DST=DST, optimals=NA, minimals=NA, HIV="nonHIV")

#make nay parameter edits here
values$months_n <- 6
# values$relapserate <- 1
# values$ltfurate_n <- values$ltfurate_sr <- 0.1


parset <- create.pars(setup = novelsetup, values = values, T, T, T)

drstate <- drout[drrow,drsetup$statenames]
newstate <- numeric(length(novelsetup$statenames)); names(newstate) <- novelsetup$statenames; newstate[drsetup$statenames] <- unlist(drstate)
s_cr <- values$cres[1]; r_cr <- values$cres[2]
novelstate <- newstate
for (name in novelsetup$statenames) if (length(grep("^S", name))==0 & length(grep("^C", name))==0 ) {  
  if (length(grep("+c+",name))==1 )   {
    if (length(grep("+.Rr+", name))==1 ) 
    {    novelstate[name] <- r_cr * newstate[str_replace(string = name, pattern = "c", replacement = "")]
    } else novelstate[name] <- s_cr * max(newstate[str_replace(string = name, pattern = "c", replacement = "")], newstate[str_replace(string = name, pattern = "c", replacement = "0")], na.rm = TRUE)
  } else    {
    if (length(grep("+.Rr+", name))==1 )
    { novelstate[name] <- (1-r_cr) * newstate[name]
    } else novelstate[name] <- (1-s_cr) * newstate[name]
  }}
    
outset6r10 <- ode(y=unlist(novelstate), times=0:10, func=dxdt, parms=parset$fullpars, do.tally=TRUE, method="adams")

var <- "tbdeaths"
plot(0:10, outset1[,var])
points(0:10, outset2[,var], pch=2)
points(0:10, outset3[,var], pch=3)
points(0:10, outset4[,var], pch=4)
points(0:10, outset6[,var], pch=6)
points(0:10, outset1r10[,var], pch=1, col='darkgreen')
points(0:10, outset2r10[,var], pch=2, col='darkgreen')
points(0:10, outset4r10[,var], pch=4, col='darkgreen')
points(0:10, outset6r10[,var], pch=6, col='darkgreen')
legend("topright",legend=c(1:4,6), pch=c(1:4,6))
points(0:10, outset1fast[,var], pch=1, col='red')
points(0:10, outset2fast[,var], pch=2, col='red')
points(0:10, outset1ltfu[,var], pch=1, col='blue')
points(0:10, outset2ltfu[,var], pch=2, col='blue')
points(0:10, outset3ltfu[,var], pch=3, col='blue')
points(0:10, outset4ltfu[,var], pch=4, col='blue')
points(0:10, outset6ltfu[,var], pch=6, col='blue')

# so the pattern (2 worse than 1) is there for ltfu 10%/mo but not 20%/mo, so it's a competition between the worse outcomes in those lost to follow up, and some other process.
# and it's even worse for a high relapse rate, so it's not the fact that they're being held inactive and nonsusceptible. 
#relapses keep dropping fom 1-4 mo, go back up for 6 mo, and are much greater for a one-mo regimen than for any others. Or if ltfu is high, relapses still drop a lot from 1 to 2 mo, but by 3 mo are the same as 1 mo regimen
#relapsers R.R0.Hn follow same pattern as relapses, they should.

# The two relevant effects I can think of are the effects that relapse happens sooner and that susceptible to reinfection happens sooner if cured. 
# I could eliminate those by having treatment last 6 months no matter what, but assigning cure vs relapse as a function of the regimen duration. 

# After compiling the two edited functions at bottom, redo:

values$months_n <- 6
parset <- create.pars(setup = novelsetup, values = values, T, T, T)
outset6fixed <- ode(y=unlist(novelstate), times=0:10, func=dxdt, parms=parset$fullpars, do.tally=TRUE, method="adams")

points(0:10, outset1fixed[,var], pch=1, col='magenta')
points(0:10, outset2fixed[,var], pch=2, col='magenta')
points(0:10, outset3fixed[,var], pch=3, col='magenta')
points(0:10, outset4fixed[,var], pch=4, col='magenta')
points(0:10, outset6fixed[,var], pch=6, col='magenta')

#Sure enough, now longer regimens are worse as they should be.

#I'll edit the dxdt function back to TPPmat version and make sure nothing changes. Yep, still looks good. 
# Then I'll edit makemat so that assignment to C vs R happens at duration rather than at 6mo for all regimens. Still ok; in fact a larger spread in the correct direction.
# And the relapsefracs(period) function itself gives output that looks just fine. 

makemat <- function(pars) 
{
  with(pars, {
    mat <- array(0, dim=c(length(Tnames), length(Tnames), length(Rnames), length(Rnames), length(Hnames), length(Hnames))); dimnames(mat) <- list("1T"=Tnames, "2T"=Tnames, "1R"=Rnames, "2R"=Rnames, "1H"=Hnames, "2H"=Hnames)
    for (jr in Rnames) for (jt in Tnames) {  mat[jt, jt, jr, jr,"Hn","Hp"] <- mat[jt, jt, jr, jr,"Hn","Hp"] + hivrate }
    for (jr in Rnames) 
    {
      mat["An","S",jr,"R0","Hn","Hn"] <- mat["An","S",jr,"R0","Hn","Hn"] + selfcurerate;  # note all cures go back to R0 to simplify later computations
      for (jt in c("Ap","Ti")) { mat[jt,"C",jr,"R0","Hn","Hn"] <- mat[jt,"C",jr,"R0","Hn","Hn"] + selfcurerate } 
      for (jh in Hnames) {  
        mat["Ln","An",jr, jr, jh, jh] <- mat["Ln","An",jr, jr, jh, jh] + reactrate[jh]; 
        mat["Lp","Ap",jr, jr, jh, jh] <- mat["Lp","Ap",jr, jr, jh, jh] + reactrate[jh]  #reactivation, at hiv-dependent rate:
        diag(mat[,,jr,jr,jh,jh]) <- diag(mat[,,jr,jr,jh,jh]) + mort[jh] #death
        diag(mat[c("An","Ap","Ti"),c("An","Ap","Ti"),jr,jr,jh,jh]) <- diag(mat[c("An","Ap","Ti"),c("An","Ap","Ti"),jr,jr,jh,jh]) + tbmort[jh] #death
        mat["R","Ap",jr,jr,jh,jh] <- mat["R","Ap",jr,jr,jh,jh] + relapserate #relapse
      }
    }
    # MODIFICATION: **************
    # once on effective treatment, progress through treatment to either relapse or cure. All including losses to follow up are assigned to C or R at 6 mo. Use only T1 for each reg.
    # ignore effect of resistance -- i.e. assume same outcomes for all 
    # (if dxdt, if on ineffective treatment (Ti), stay there 6 mo then restart treatment.)
    rduration <- numeric(0)
    rduration[1] <- 0.03*  0.09806163 + 0.97*(0.01903855)
    rduration[2] <-  0.03*(0.34746851  + 0.05438356) + 0.94*  0.01903855  
    rduration[3] <-  0.03*(0.56497901  + 0.09806163+  0.04260189) + 0.91*  0.01903855 
    rduration[4] <-   0.03*(0.67373426  +0.12198392 + 0.07413933+  0.03671106) + 0.88* 0.01903855
    rduration[6] <- 0.03*(0.78248950 +0.34746851  +0.11400982  +0.08211343 +2*0.04260189) + 0.82*(0.01903855)

    for (jh in 1:length(Hnames))
    {
      for (jr in Rnames) for (reg in regimens[usereg])
      {
        relapses <- ifelse(durations[reg]>0.5, relapsepoor*poor_r, max(poormat[reg,1]/poormat['n',1]*rduration[round(12*durations[reg])]))
        period <- paste0("T", reg, 1)
        mat[period,"R",jr,jr,jh,jh] <- 
            mat[period,"R",jr,jr,jh,jh] + (12*durations[reg])*relapses
        mat[period,"C",jr,"R0",jh,jh] <- 
            mat[period,"C",jr,"R0",jh,jh] + (12*durations[reg])*(1-relapses)
        }
      }
    return(mat)
  })
}
