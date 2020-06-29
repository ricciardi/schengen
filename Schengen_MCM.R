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

outcome.vars <- c("CBWbord","CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_0","seekdur_1_2","seekdur_3more")

for(o in outcome.vars){

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
  
  # Block resampling with fixed block lengths of length l)
  source("MCEstBoot.R")
  
  boot <- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X, X.hat=outcomes.cbw$X.hat, covars=FALSE, rev=TRUE, R=10, parallel = "multicore", l=bopt, sim = "fixed") 
  saveRDS(boot, paste0("results/boot-cbw-",o,".rds"))
  
  t0 <- which(colnames(outcomes.cbw$M)=="20091")
  treat_indices_order <- outcomes.lm$treated
  
  boot.trajectory <- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X, X.hat=outcomes.cbw$X.hat, treat_indices_order=treat_indices_order, covars=FALSE, rev=TRUE, t0=t0, R=1000, parallel = "multicore", l=bopt, sim = "fixed") 
  saveRDS(boot.trajectory, paste0("results/boot-trajectory-cbw-",o,".rds")) # bootstrap for ATT trajectory
  
  print(paste0("Analysis 1, outcome",o, "ATT:",boot.trajectory$t0,"CI:",boot.ci(boot.trajectory)))
  
  ## Analysis 2: ST vs NT (forward, X=LM)
  
  if(o %in% c("CBWbord","CBWbordEMPL")) next
  
  print(paste0("Estimates for Analysis 2, outcome:",o))
  
  outcomes.lm <- readRDS(paste0("data/outcomes-lm-",o,".rds"))
  
  # Get treatment effect estimates
  
  mc.estimates.lm <- MCEst(outcomes.lm, rev=FALSE, covars=FALSE) 
  saveRDS(mc.estimates.lm, paste0("results/mc-estimates-lm-",o,".rds"))
  
  # Get optimal stationary bootstrap lengths

  bopt <- b.star(t(outcomes.lm$M),round=TRUE)[,1]
  
  # Block resampling with fixed block lengths of length l)
  
  boot <- tsboot(tseries=ts(t(outcomes.lm$M)), MCEstBoot, mask=outcomes.lm$mask, W=outcomes.lm$W, X=outcomes.lm$X,X.hat=outcomes.lm$X.hat,covars=FALSE,rev=FALSE,R=1000, parallel = "multicore", l=bopt, sim = "fixed") 
  saveRDS(boot, paste0("results/boot-lm-",o,".rds"))
  
  # Get p-values

  t0 <- which(colnames(outcomes.lm$M)=="20091")
  
  treat_indices_order <- outcomes.lm$treated
  
  boot.trajectory <- tsboot(tseries=ts(t(outcomes.lm$M)), MCEstBoot, mask=outcomes.lm$mask, W=outcomes.lm$W, X=outcomes.lm$X, X.hat=outcomes.lm$X.hat, treat_indices_order=treat_indices_order, covars=FALSE, rev=FALSE, t0=t0, R=1000, parallel = "multicore", l=bopt, sim = "fixed") 
  saveRDS(boot.trajectory, paste0("results/boot-trajectory-lm-",o,".rds")) # bootstrap for ATT trajectory
  
  print(paste0("Analysis 1, outcome",o, "ATT:",boot.trajectory$t0,"CI:",boot.ci(boot.trajectory)))
}