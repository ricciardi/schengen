###################################
# Placebo MC estimates with covariates #
###################################

## Loading Source files
library(MCPanel)
library(glmnet)
library(ggplot2)

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
  
  # Eastern cluster
  print(paste0("Estimates for Analysis 1, outcome:",o))
  
  outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  # Discard post-treatment periods
  outcomes.cbw.placebo <- outcomes.cbw
  outcomes.cbw.placebo$M <- outcomes.cbw$M[,1:which(colnames(outcomes.cbw$M)=="20091")]
  outcomes.cbw.placebo$W <- outcomes.cbw$W[,1:which(colnames(outcomes.cbw$W)=="20091")]
  outcomes.cbw.placebo$mask <- outcomes.cbw$mask[,1:which(colnames(outcomes.cbw$mask)=="20091")]
  
  # Get p-values
  source("MCEst.R")
  source("ChernoTest.R")
  
  t_final_placebo <- ncol(outcomes.cbw.placebo$M ) # all periods 
  
  taus <- 1:length((4:t_final_placebo))
  
  treat_indices_order <- outcomes.cbw.placebo$treated
  
  moving.block.placebo <- lapply(taus, function(t){
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.cbw.placebo[c("M","mask","W")], ns=1000, treat_indices_order=treat_indices_order, permtype="moving.block",t0=t0_placebo,rev=TRUE,covars=TRUE)
    })
  saveRDS(moving.block.placebo,paste0("results/moving-block-placebo-cbw-",o,"-covars.rds"))
  
  iid.block.placebo <- lapply(taus, function(t){
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.cbw.placebo[c("M","mask","W")], ns=1000, treat_indices_order=treat_indices_order, permtype="iid.block",t0=t0_placebo,rev=TRUE,covars=TRUE)
    })
  saveRDS(iid.block.placebo,paste0("results/iid-block-placebo-cbw-",o,"-covars.rds"))
  
  iid.placebo <- lapply(taus, function(t){
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.cbw.placebo[c("M","mask","W")],ns=1000, treat_indices_order=treat_indices_order, permtype="iid",t0=t0_placebo,rev=TRUE,covars=TRUE)
    })
  saveRDS(iid.placebo,paste0("results/iid-placebo-cbw-",o,"-covars.rds"))
  
  ## Analysis 2: ST vs NT (forward, X=LM)
  
  # Eastern cluster
  print(paste0("Estimates for Analysis 2, outcome:",o))
  
  outcomes.lm <- readRDS(paste0("data/outcomes-lm-",o,".rds"))
  
  # Discard post-treatment periods
  outcomes.lm.placebo <- outcomes.lm
  outcomes.lm.placebo$M <- outcomes.lm$M[,1:which(colnames(outcomes.lm$M)=="20091")-1]
  outcomes.lm.placebo$W <- outcomes.lm$W[,1:which(colnames(outcomes.lm$W)=="20091")-1]
  outcomes.lm.placebo$mask <- outcomes.lm$mask[,1:which(colnames(outcomes.lm$mask)=="20091")-1]
  
  # Get p-values
  
  t_final_placebo <- ncol(outcomes.lm.placebo$M ) # all periods 
  
  taus <- 1:length((4:t_final_placebo))
  
  treat_indices_order <- outcomes.lm.placebo$treated
  
  moving.block.placebo <- lapply(taus, function(t){
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.lm.placebo[c("M","mask","W")], ns=1000, treat_indices_order=treat_indices_order, permtype="moving.block",t0=t0_placebo,rev=FALSE,covars=TRUE)
    })
  saveRDS(moving.block.placebo,paste0("results/moving-block-placebo-lm-",o,"-covars.rds"))
  
  iid.block.placebo <- lapply(taus, function(t){
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.lm.placebo[c("M","mask","W")], ns=1000, treat_indices_order=treat_indices_order, permtype="iid.block",t0=t0_placebo,rev=FALSE,covars=TRUE)
    })
  saveRDS(iid.block.placebo,paste0("results/iid-block-placebo-lm-",o,"-covars.rds"))
  
  iid.placebo <- lapply(taus, function(t){
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.lm.placebo[c("M","mask","W")],ns=1000, treat_indices_order=treat_indices_order, permtype="iid",t0=t0_placebo,rev=FALSE,covars=TRUE)
    })
  saveRDS(iid.placebo,paste0("results/iid-placebo-lm-",o,"-covars.rds"))
}