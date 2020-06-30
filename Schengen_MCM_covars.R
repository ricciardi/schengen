###################################
#  MC estimates with covariates#
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
  mc.estimates.cbw <- MCEst(outcomes.cbw, rev=TRUE, covars=TRUE) 
  saveRDS(mc.estimates.cbw, paste0("results/mc-estimates-cbw-",o,"-covars.rds"))
  
  # Get optimal stationary bootstrap lengths
  source("PolitisWhite.R")
  
  bopt <- b.star(t(outcomes.cbw$M),round=TRUE)[,1]
  
  # Block resampling with fixed block lengths of length l)
  source("MCEstBoot.R")
  
  boot.cbw <- tsboot(tseries=t(outcomes.cbw$M), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X,X.hat=outcomes.cbw$X.hat, covars=TRUE,R=1000, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot.cbw, paste0("results/boot-cbw-",o,"-covars.rds"))
  
  # Bootstrap for ATT trajectory 
  
  t0_eastern.cbw <- which(colnames(outcomes.cbw$M)=="20111") # earliest combined treatment
  t0_swiss.cbw <- which(colnames(outcomes.cbw$M)=="20091") # earliest combined treatment

  boot.trajectory.eastern <- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X,X.hat=outcomes.cbw$X.hat,eastern=outcomes.cbw$eastern, covars=TRUE, rev=TRUE, t0=t0_eastern.cbw, R=1000, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot.trajectory.eastern, paste0("results/boot-trajectory-eastern-cbw-",o,"-covars.rds")) 
  
  print(paste0("Analysis 1, Eastern cluster, Outcome",o))
  print(mean(boot.trajectory.eastern$t0)) # test statistic S
  print(boot.ci(boot.trajectory.eastern, type=c("basic","norm")))
  
  boot.trajectory.swiss <- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X,X.hat=outcomes.cbw$X.hat,swiss=outcomes.cbw$swiss, covars=TRUE, rev=TRUE, t0=t0_swiss.cbw, R=1000, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot.trajectory.swiss, paste0("results/boot-trajectory-swiss-cbw-",o,"-covars.rds")) 
  
  print(paste0("Analysis 1, swiss cluster, Outcome",o))
  print(mean(boot.trajectory.swiss$t0)) # post-period att
  print(boot.ci(boot.trajectory.swiss, type=c("basic","norm")))

  ## Analysis 2: ST vs NT (forward, X=LM)
  
  if(o %in% c("N_CBWbord","CBWbord","CBWbordEMPL")) next
  
  print(paste0("Estimates for Analysis 2, outcome:",o))
  
  outcomes.lm <- readRDS(paste0("data/outcomes-lm-",o,".rds"))
  
  # Get treatment effect estimates
  
  mc.estimates.lm <- MCEst(outcomes.lm, rev=FALSE, covars=TRUE) 
  saveRDS(mc.estimates.lm, paste0("results/mc-estimates-lm-",o,"-covars.rds"))
  
  # Get optimal stationary bootstrap lengths

  bopt <- b.star(t(outcomes.lm$M),round=TRUE)[,1]
  
  # Block resampling with fixed block lengths of length l)
  
  boot.lm <- tsboot(tseries=ts(t(outcomes.lm$M)), MCEstBoot, mask=outcomes.lm$mask, W=outcomes.lm$W, X=outcomes.lm$X,X.hat=outcomes.lm$X.hat, covars=TRUE,rev=FALSE, R=1000, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot.lm, paste0("results/boot-lm-",o,"-covars.rds"))
  
  # Get CIs for trajectory
  
  t0_eastern.lm <- which(colnames(outcomes.cbw$M)=="20081") # schengen start
  t0_swiss.lm <- which(colnames(outcomes.cbw$M)=="20072") # FoM start

  boot.trajectory.eastern <- tsboot(tseries=ts(t(outcomes.lm$M)), MCEstBoot, mask=outcomes.lm$mask, W=outcomes.lm$W, eastern=outcomes.lm$eastern,X=outcomes.lm$X,X.hat=outcomes.lm$X.hat, covars=TRUE, rev=FALSE, t0=t0_eastern.lm, R=1000, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot.trajectory.eastern, paste0("results/boot-trajectory-eastern-lm-",o,"-covars.rds")) 
  
  print(paste0("Analysis 1, Eastern cluster, Outcome",o))
  print(mean(boot.trajectory.eastern$t0)) # test statistic S
  print(boot.ci(boot.trajectory.eastern, type=c("basic","norm")))
  
  boot.trajectory.swiss <- tsboot(tseries=ts(t(outcomes.lm$M)), MCEstBoot, mask=outcomes.lm$mask, W=outcomes.lm$W, swiss=outcomes.lm$swiss, X=outcomes.lm$X,X.hat=outcomes.lm$X.hat,covars=TRUE, rev=FALSE, t0=t0_swiss.lm, R=1000, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot.trajectory.swiss, paste0("results/boot-trajectory-swiss-lm-",o,"-covars.rds")) 
  
  print(paste0("Analysis 1, swiss cluster, Outcome",o))
  print(mean(boot.trajectory.swiss$t0)) # post-period att
  print(boot.ci(boot.trajectory.swiss, type=c("basic","norm")))
}