###################################
# Placebo MC estimates #
###################################

## Loading Source files
library(MCPanel)
library(glmnet)
library(ggplot2)
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

  ## Analysis 1: ST vs AT (retrospective, X=CBW) 
  
  print(paste0("Estimates for Analysis 1, outcome:",o))
  
  outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  # Use pre-treatment for LT (no missing values)
  outcomes.cbw.placebo <- outcomes.cbw
  outcomes.cbw.placebo$mask <- outcomes.cbw$mask[,(which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask))]  # all zeros
  outcomes.cbw.placebo$M <- outcomes.cbw$M[,(which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask))]
  outcomes.cbw.placebo$W <- outcomes.cbw$W[,(which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask))]
  
  # Get p-values
  source("MCEst.R")
  source("ChernoTest.R")
  
  t_final_placebo <- ncol(outcomes.cbw.placebo$M ) # all periods
  print(t_final_placebo) 
  
  t0_placebo <- c(round(t_final_placebo/10), round(t_final_placebo/8), round(t_final_placebo/5), round(t_final_placebo/3)) # n pre-treatment periods
  
  iid.block.placebo.eastern <- lapply(t0_placebo, function(t){
    ChernoTest(outcomes=outcomes.cbw.placebo[c("M","mask","W","X","X.hat")], ns=1000, q=1,treat_indices_order=outcomes.cbw$eastern, permtype="iid.block",t0=t,rev=TRUE,covars=FALSE)
    
  })
  saveRDS(iid.block.placebo.eastern,paste0("results/iid-block-placebo-cbw-eastern",o,".rds"))

    iid.block.placebo.swiss <- lapply(t0_placebo, function(t){
    ChernoTest(outcomes=outcomes.cbw.placebo[c("M","mask","W","X","X.hat")], ns=1000, q=1,treat_indices_order=outcomes.cbw$swiss, permtype="iid.block",t0=t,rev=TRUE,covars=FALSE)
    
  })
  saveRDS(iid.block.placebo.swiss,paste0("results/iid-block-placebo-cbw-swiss",o,".rds"))
}