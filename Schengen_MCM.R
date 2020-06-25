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
  print(o)
  
  ## Analysis 1: ST vs AT (retrospective, X=CBW) 
  
  print(paste0("Estimates for Analysis 1, outcome:",o))
  
  outcomes.cbw <- readRDS(paste0("data/outcomes-cbw",o,".rds"))
  
  # Get treatment effect estimates
    
  source('MCEst.R')
  mc.estimates.cbw <- MCEst(outcomes.cbw, covars=FALSE) 
  saveRDS(mc.estimates.cbw, paste0("results/mc-estimates-cbw",o,".rds"))
  
  # Get optimal stationary bootstrap lengths
  source("PolitisWhite.R")
  
  bopt <- b.star(t(outcomes.cbw$M),round=TRUE)[,1]
  
  # Block resampling with fixed block lengths of length l)
  source("MCEstBoot.R")
  
  boot <- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, treated=outcomes.cbw$treated, control=outcomes.cbw$control, covars=FALSE,R=1000, parallel = "multicore", l=bopt, sim = "fixed") 
  saveRDS(boot, paste0("results/boot-cbw",o,".rds"))
  
  # Get p-values
  source("ChernoTest.R")
  
  t0 <- which(colnames(outcomes.cbw$M)=="20091")
  
  treat_indices_order <- outcomes.cbw$treated
  
  iid.block <- ChernoTest(outcomes=outcomes.cbw[c("M","mask","W")], ns=1000, treat_indices_order=treat_indices_order, permtype="iid.block",t0=t0,rev=TRUE,covars=FALSE)
  saveRDS(iid.block,paste0("results/iid-block-cbw",o,".rds"))
  
  iid <- ChernoTest(outcomes=outcomes.cbw[c("M","mask","W")], ns=1000, treat_indices_order=treat_indices_order, permtype="iid",t0=t0,rev=TRUE,covars=FALSE)
  saveRDS(iid,paste0("results/iid-cbw",o,".rds"))
  
  moving.block <- ChernoTest(outcomes=outcomes.cbw[c("M","mask","W")], ns=1000, treat_indices_order=treat_indices_order, permtype="moving.block",t0=t0,rev=TRUE,covars=FALSE)
  saveRDS(moving.block,paste0("results/moving-block-cbw",o,".rds"))

  ## Analysis 2: ST vs NT (forward, X=LM)
  
  print(paste0("Estimates for Analysis 2, outcome:",o))
  
  outcomes.lm <- readRDS(paste0("data/outcomes-lm",o,".rds"))
  
  # Get treatment effect estimates
  
  mc.estimates.lm <- MCEst(outcomes.lm, covars=FALSE) 
  saveRDS(mc.estimates.lm, paste0("results/mc-estimates-lm",o,".rds"))
  
  # Get optimal stationary bootstrap lengths

  bopt <- b.star(t(outcomes.lm$M),round=TRUE)[,1]
  
  # Block resampling with fixed block lengths of length l)
  
  boot <- tsboot(tseries=ts(t(outcomes.lm$M)), MCEstBoot, mask=outcomes.lm$mask, W=outcomes.lm$W, treated=outcomes.lm$treated, control=outcomes.lm$control, covars=FALSE,R=1000, parallel = "multicore", l=bopt, sim = "fixed") 
  saveRDS(boot, paste0("results/boot-lm",o,".rds"))
  
  # Get p-values

  t0 <- which(colnames(outcomes.lm$M)=="20091")
  
  treat_indices_order <- outcomes.lm$treated
  
  iid.block <- ChernoTest(outcomes=outcomes.lm[c("M","mask","W")], ns=1000, treat_indices_order=treat_indices_order, permtype="iid.block",t0=t0,covars=FALSE)
  saveRDS(iid.block,paste0("results/iid-block-lm",o,".rds"))
  
  moving.block <- ChernoTest(outcomes=outcomes.lm[c("M","mask","W")], ns=1000, treat_indices_order=treat_indices_order, permtype="moving.block",t0=t0,covars=FALSE)
  saveRDS(moving.block,paste0("results/moving-block-lm",o,".rds"))
  
  iid <- ChernoTest(outcomes=outcomes.lm[c("M","mask","W")], ns=1000, treat_indices_order=treat_indices_order, permtype="iid",t0=t0,covars=FALSE)
  saveRDS(iid,paste0("results/iid-lm",o,".rds"))
}