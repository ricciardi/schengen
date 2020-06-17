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
  mc.estimates.cbw.eastern <- MCEst(outcomes.cbw.eastern, imputed=FALSE, covars=NULL) # always-treated and switch-treated
  saveRDS(mc.estimates.cbw.eastern, "results/mc-estimates-cbw-eastern.rds")
  
  # Get optimal stationary bootstrap lengths
  source("PolitisWhite.R")
  
  bopt <- b.star(t(outcomes.cbw.eastern$M),round=TRUE)[,1]
  
  # Block resampling with fixed block lengths of length l)
  source("MCEstBoot.R")
  
  boot <- tsboot(tseries=ts(t(outcomes.cbw.eastern$M)), MCEstBoot, M.missing=outcomes.cbw.eastern$M.missing, mask=outcomes.cbw.eastern$mask, imputed=TRUE,covars=NULL,R=1000, parallel = "multicore", l=bopt, sim = "fixed") 
  saveRDS(boot, "results/boot.rds")
  
  pvals <- FALSE
  if(pvals){
    # Get p-values
    source("ChernoTest.R")
    
    t0 <- which(colnames(outcomes.cbw.eastern$M)=="20081")
    
    treat_indices_order <- outcomes.cbw.eastern$treated
    
    moving.block <- ChernoTest(outcomes=outcomes.cbw.eastern[c("M","M.missing","mask")], ns=500, treat_indices_order=treat_indices_order, permtype="moving.block",t0=t0,imputed=TRUE,covars=NULL)
    saveRDS(moving.block,"results/moving_block.rds")
    
    iid.block <- ChernoTest(outcomes=outcomes.cbw.eastern[c("M","M.missing","mask")], ns=500, treat_indices_order=treat_indices_order, permtype="iid.block",t0=t0,imputed=TRUE,covars=NULL)
    saveRDS(iid.block,"results/iid_block.rds")
    
    iid <- ChernoTest(outcomes=outcomes.cbw.eastern[c("M","M.missing","mask")],ns=500, treat_indices_order=treat_indices_order, permtype="iid",t0=t0,imputed=TRUE,covars=NULL)
    saveRDS(iid,"results/iid.rds")
    
    mean(iid$real_att)
    iid$p
    iid.block$p
    moving.block$p
  }
}