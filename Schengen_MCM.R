###################################
#  MC estimates #
###################################

## Libraries
library(MCPanel)
library(boot)
library(glmnet)

# Setup parallel processing 
library(parallel)
library(doParallel)

cores <- detectCores()

cl <- parallel::makeForkCluster(cores)

doParallel::registerDoParallel(cores) # register cores (<p)

RNGkind("L'Ecuyer-CMRG") # ensure random number generation

outcome.vars <- c("N_CBWbord","CBWbord","CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_0","seekdur_1_2","seekdur_3more")

for(o in outcome.vars){
  print(o)

  ## Analysis 1: ST vs AT (retrospective, X=CBW) 
  
  print(paste0("Estimates for Analysis 1, outcome:",o))
  
  outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  # Get treatment effect estimates
    
  source('MCEst.R')
  mc.estimates.cbw <- MCEst(outcomes.cbw, rev=TRUE, covars=FALSE) 
  saveRDS(mc.estimates.cbw, paste0("results/mc-estimates-cbw-",o,".rds"))
  
  # Get optimal stationary bootstrap lengths
  source("PolitisWhite.R")
  
  bopt <- b.star(t(outcomes.cbw$M),round=TRUE)[,1]
  
  # Bootstrap for ATT trajectory 
  source("MCEstBoot.R")
  
  t0_eastern.cbw <- which(colnames(outcomes.cbw$M)=="20111") # earliest combined treatment
  t0_swiss.cbw <- which(colnames(outcomes.cbw$M)=="20091") # earliest combined treatment

  boot.trajectory.eastern <- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, eastern=outcomes.cbw$eastern, covars=FALSE, rev=TRUE, t0=t0_eastern.cbw, hyp=FALSE, R=999, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot.trajectory.eastern, paste0("results/boot-trajectory-eastern-cbw-",o,".rds")) 
  
  boot.trajectory.swiss <- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, swiss=outcomes.cbw$swiss, covars=FALSE, rev=TRUE, t0=t0_swiss.cbw, hyp=FALSE, R=999, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot.trajectory.swiss, paste0("results/boot-trajectory-swiss-cbw-",o,".rds")) 
}