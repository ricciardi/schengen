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

cores <- 2#detectCores()

cl <- parallel::makeForkCluster(cores)

doParallel::registerDoParallel(cores) # register cores (<p)

RNGkind("L'Ecuyer-CMRG") # ensure random number generation

outcomes <- c("CBWbord","CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_0","seekdur_1_2","seekdur_3more")

for(o in outcomes){
  print(o)
  
  ## Analysis 1: ST vs AT (retrospective, X=CBW) 
  
  # Eastern cluster
  
  outcomes.cbw.eastern <- readRDS(paste0("data/outcomes-cbw-eastern-",o,".rds"))
  
  outcomes.cbw.eastern$mask[rownames(outcomes.cbw.eastern$mask)%in%outcomes.cbw.eastern$treated,] <- abs(outcomes.cbw.eastern$mask[rownames(outcomes.cbw.eastern$mask)%in%outcomes.cbw.eastern$treated,]-1) # estimate Y(1)_LT,pre
  
  # Get treatment effect estimates
    
  source('MCEst.R')
  mc.estimates.cbw.eastern <- MCEst(outcomes.cbw.eastern, covars=NULL) 
  saveRDS(mc.estimates.cbw.eastern, "results/mc-estimates-cbw-eastern.rds")
  
  # Get optimal stationary bootstrap lengths
  source("PolitisWhite.R")
  
  bopt <- b.star(t(outcomes.cbw.eastern$M),round=TRUE)[,1]
  
  # Block resampling with fixed block lengths of length l)
  source("MCEstBoot.R")
  
  boot <- tsboot(tseries=ts(t(outcomes.cbw.eastern$M)), MCEstBoot, mask=outcomes.cbw.eastern$mask, W=outcomes.cbw.eastern$W, covars=NULL,R=1000, parallel = "multicore", l=bopt, sim = "fixed") 
  saveRDS(boot, "results/boot-cbw-eastern.rds")
  
  # Get p-values
  source("ChernoTest.R")
  
  t0 <- which(colnames(outcomes.cbw.eastern$M)=="20081")
  
  treat_indices_order <- outcomes.cbw.eastern$treated
  
  iid.block <- ChernoTest(outcomes=outcomes.cbw.eastern[c("M","mask","W")], ns=1000, treat_indices_order=treat_indices_order, permtype="iid.block",t0=t0,covars=NULL)
  saveRDS(iid.block,"results/iid-block-cbw-eastern.rds")
  
  # Swiss cluster
  
  outcomes.cbw.swiss <- readRDS(paste0("data/outcomes-cbw-swiss-",o,".rds"))
  
  outcomes.cbw.swiss$mask[rownames(outcomes.cbw.swiss$mask)%in%outcomes.cbw.swiss$treated,] <- abs(outcomes.cbw.swiss$mask[rownames(outcomes.cbw.swiss$mask)%in%outcomes.cbw.swiss$treated,]-1) # estimate Y(1)_LT,pre
  
  # Get treatment effect estimates
  
  mc.estimates.cbw.swiss <- MCEst(outcomes.cbw.swiss, covars=NULL)
  saveRDS(mc.estimates.cbw.swiss, "results/mc-estimates-cbw-swiss.rds")
  
  # Get optimal stationary bootstrap lengths
  source("PolitisWhite.R")
  
  bopt <- b.star(t(outcomes.cbw.swiss$M),round=TRUE)[,1]
  
  # Block resampling with fixed block lengths of length l)
  
  boot <- tsboot(tseries=ts(t(outcomes.cbw.swiss$M)), MCEstBoot, mask=outcomes.cbw.swiss$mask, W=outcomes.cbw.swiss$W, covars=NULL,R=1000, parallel = "multicore", l=bopt, sim = "fixed") 
  saveRDS(boot, "results/boot-cbw-swiss.rds")
  
  # Get p-values
  
  t0 <- which(colnames(outcomes.cbw.swiss$M)=="20072")
  
  treat_indices_order <- outcomes.cbw.swiss$treated
  
  iid.block <- ChernoTest(outcomes=outcomes.cbw.swiss[c("M","mask","W")], ns=1000, treat_indices_order=treat_indices_order, permtype="iid.block",t0=t0,covars=NULL)
  saveRDS(iid.block,"results/iid-block-cbw-swiss.rds")

  ## Analysis 2: ST vs NT (forward, X=LM)
  
  # Eastern cluster
  
  outcomes.lm.eastern <- readRDS(paste0("data/outcomes-lm-eastern-",o,".rds"))
  
  # Get treatment effect estimates
  
  mc.estimates.lm.eastern <- MCEst(outcomes.lm.eastern, covars=NULL) 
  saveRDS(mc.estimates.lm.eastern, "results/mc-estimates-lm-eastern.rds")
  
  # Get optimal stationary bootstrap lengths

  bopt <- b.star(t(outcomes.lm.eastern$M),round=TRUE)[,1]
  
  # Block resampling with fixed block lengths of length l)
  
  boot <- tsboot(tseries=ts(t(outcomes.lm.eastern$M)), MCEstBoot, mask=outcomes.lm.eastern$mask, W=outcomes.lm.eastern$W, covars=NULL,R=1000, parallel = "multicore", l=bopt, sim = "fixed") 
  saveRDS(boot, "results/boot-lm-eastern.rds")
  
  # Get p-values

  t0 <- which(colnames(outcomes.lm.eastern$M)=="20081")
  
  treat_indices_order <- outcomes.lm.eastern$treated
  
  iid.block <- ChernoTest(outcomes=outcomes.lm.eastern[c("M","mask","W")], ns=1000, treat_indices_order=treat_indices_order, permtype="iid.block",t0=t0,covars=NULL)
  saveRDS(iid.block,"results/iid-block-lm-eastern.rds")
  
  # Swiss cluster
  
  outcomes.lm.swiss <- readRDS(paste0("data/outcomes-lm-swiss-",o,".rds"))
  
  # Get treatment effect estimates
  
  mc.estimates.lm.swiss <- MCEst(outcomes.lm.swiss, covars=NULL) 
  saveRDS(mc.estimates.lm.swiss, "results/mc-estimates-lm-swiss.rds")
  
  # Get optimal stationary bootstrap lengths
  
  bopt <- b.star(t(outcomes.lm.swiss$M),round=TRUE)[,1]
  
  # Block resampling with fixed block lengths of length l)
  
  boot <- tsboot(tseries=ts(t(outcomes.lm.swiss$M)), MCEstBoot, mask=outcomes.lm.swiss$mask, W=outcomes.lm.swiss$W, covars=NULL,R=1000, parallel = "multicore", l=bopt, sim = "fixed") 
  saveRDS(boot, "results/boot-lm-swiss.rds")
  
  # Get p-values
  
  t0 <- which(colnames(outcomes.lm.swiss$M)=="20072")
  
  treat_indices_order <- outcomes.lm.swiss$treated
  
  iid.block <- ChernoTest(outcomes=outcomes.lm.swiss[c("M","mask","W")], ns=1000, treat_indices_order=treat_indices_order, permtype="iid.block",t0=t0,covars=NULL)
  saveRDS(iid.block,"results/iid-block-lm-swiss.rds")
}