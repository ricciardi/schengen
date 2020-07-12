###################################
#  MC estimates with covariates#
###################################

## Libraries
library(MCPanel)
library(boot)

# Setup parallel processing 
library(parallel)
library(doParallel)

cores <- detectCores()

cl <- parallel::makeForkCluster(cores)

doParallel::registerDoParallel(cores) # register cores (<p)

RNGkind("L'Ecuyer-CMRG") # ensure random number generation

outcome.vars <- c("CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_0","seekdur_1_2","seekdur_3more")

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
  
  # Bootstrap for per-period effects
  source("MCEstBoot.R")
  
  boot <- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X,X.hat=outcomes.cbw$X.hat,covars=TRUE, rev=TRUE, R=999, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot, paste0("results/boot-cbw-",o,"-covars.rds")) 
  
  # Bootstrap by group

  boot.eastern<- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X,X.hat=outcomes.cbw$X.hat,eastern=outcomes.cbw$eastern, covars=TRUE, rev=TRUE, R=999, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot.eastern, paste0("results/boot-eastern-cbw-",o,"-covars.rds")) 
  
  boot.swiss<- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X,X.hat=outcomes.cbw$X.hat,swiss=outcomes.cbw$swiss, covars=TRUE, rev=TRUE, R=999, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot.swiss, paste0("results/boot-swiss-cbw-",o,"-covars.rds")) 
  
  boot.control<- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X,X.hat=outcomes.cbw$X.hat,control=outcomes.cbw$control, covars=TRUE, rev=TRUE, R=999, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot.control, paste0("results/boot-control-cbw-",o,"-covars.rds")) 
}