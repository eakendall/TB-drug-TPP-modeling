# linear declind in parameter value (beta and reactivation, separately) to get a 2% and a 5% annual decline in incidence
# then, vary all params in the "neg direction" simultaneously from ~-25 to +5%, OR,
# go through all the parameters and make sure reasonable changes in those don't change TB incidence by much
# , or just look at PRCCs for incidence to justify the choice of the first two


# first look at PRCCs, for incidence projection:
library("sensitivity")
prccinc <- pcc(X = novelwide1[,41:75], y= alldrDST[,"inc10"], rank=TRUE ) 
cbind(rownames(prccinc$PRCC), prccinc$PRCC)[rev(order(abs(prccinc$PRCC))),]
# most influential are self cure rate, diagnostic rate, relapse vs failure, tb and baseline mortality, MDR acquisition and transmissibility, reactivation rate, 

# but wait, need to look at uncalibrated sensitivity. can pull from optimats? Or can do one-way sensitivity analyses? 

# set up one-ways: 
midvalues <- set.values(pessimistic = TRUE)
midvalues$cal$beta <- 4.5; midvalues$cal$hivrate <- 0.0004
midparams <- create.pars(values = midvalues, pessimistic = TRUE, DRera = FALSE, treatSL = FALSE, treatnovel = FALSE)
e <- equilib(pars = midparams, tol=0.5) # run to equilib 
e$log[nrow(e$log),"inc"] 
e$log[nrow(e$log),"prev"] 
e$log[nrow(e$log),"coprev"]
state <- e$log[nrow(e$log),2:length(dssetup$statenames)]

write(c("param", "index", "lowinc", "highinc"), file = "incsensis.csv", ncolumns = 4, append = FALSE, sep=",")

j <- 0
for ( m in 3:4 )
for ( n in 1:length(midvalues[[m]]) )
for ( j in 1:length(unlist(midvalues[[m]][n])) )
{
  lowvalues <- highvalues <- midvalues
  lowvalues[[m]][[n]][j] <- midvalues[[m]][[n]][j] * 0.75
  highvalues[[m]][[n]][j] <- midvalues[[m]][[n]][j] * 1.25
  lowparams <- create.pars(values = lowvalues, pessimistic = TRUE, DRera = FALSE, treatSL = FALSE, treatnovel = FALSE)
  highparams <- create.pars(values = highvalues, pessimistic = TRUE, DRera = FALSE, treatSL = FALSE, treatnovel = FALSE)

  elow <- equilib(pars = lowparams, tol=0.5) 
  ehigh <- equilib(pars = highparams, tol=0.5) 

  write(c(names(midparams$values[[m]])[n], j, elow$log[nrow(elow$log),"inc"], ehigh$log[nrow(ehigh$log),"inc"]), file = "incsensis.csv", ncolumns = 4, append = TRUE, sep=",")
}

r <- read.csv("incsensis.csv", header=TRUE)
ratios <- (r[,"highinc"]-r[,"lowinc"])/e$log[nrow(e$log),"inc"] 
names(ratios) <- r[,"param"]
ratios
-0.05*0.5/ratios[order(abs(ratios))] # change needed for 5% decline
# choose those that make largest impact: reactrate1, dxrate1, rapidprog1, and beta

# I edited dxdt to vary a parameter if ode is called with argument paramvary=list(paramname, annualchange)
# For params I want to vary, I'll need to figure out annual reduction (aiming for 2-5%/year, given that ~50% change gives r): 
beta -0.02, rapidprog -0.02, dxrate +0.05, reactrate -0.04
# Then, for a subset of DR runs (say, start with idr==1), I'll redo allminopt, allbut, and only , with paramvary in the ode calls

r <- read.csv("Allminopt.reactrate-0.02_DSDSTall_India_20160419p.1.csv")
r <- read.csv("Only.rapidprog-0.015_DSDSTall_India_20160419p.1.csv")
summary(r$inc0allmin)
summary(r$inc10allmin)
125/147
