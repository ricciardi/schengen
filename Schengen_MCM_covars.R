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
  
  # Bootstrap for ATT trajectory CIs
  source("MCEstBoot.R")
  
  t0_eastern.cbw <- which(colnames(outcomes.cbw$M)=="20111") # earliest combined treatment
  t0_swiss.cbw <- which(colnames(outcomes.cbw$M)=="20091") # earliest combined treatment

  boot.trajectory.eastern <- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X,X.hat=outcomes.cbw$X.hat,eastern=outcomes.cbw$eastern, covars=TRUE, rev=TRUE, t0=t0_eastern.cbw, R=999, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot.trajectory.eastern, paste0("results/boot-trajectory-eastern-cbw-",o,"-covars.rds")) 
  
  boot.trajectory.swiss <- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X,X.hat=outcomes.cbw$X.hat,swiss=outcomes.cbw$swiss, covars=TRUE, rev=TRUE, t0=t0_swiss.cbw, R=999, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot.trajectory.swiss, paste0("results/boot-trajectory-swiss-cbw-",o,"-covars.rds")) 
  
  # Get p-values
  source("ChernoTest.R")
  
  iid.block.eastern <- ChernoTest(outcomes=outcomes.cbw[c("M","mask","W","X","X.hat")], ns=1000, q=1, t.stat=boot.trajectory.eastern$t0, treat_indices_order=outcomes.cbw$eastern, permtype="iid.block",t0=t0_eastern.cbw,rev=TRUE,covars=TRUE,bopt=round(mean(bopt)))
  saveRDS(iid.block.eastern,paste0("results/iid-block-cbw-eastern-",o,"-covars.rds"))
  
  iid.block.swiss <- ChernoTest(outcomes=outcomes.cbw[c("M","mask","W","X","X.hat")], ns=1000, q=1, t.stat=boot.trajectory.swiss$t0, treat_indices_order=outcomes.cbw$swiss, permtype="iid.block",t0=t0_swiss.cbw,rev=TRUE,covars=TRUE,bopt=round(mean(bopt)))
  saveRDS(iid.block.swiss,paste0("results/iid-block-cbw-swiss-",o,"-covars.rds"))

}