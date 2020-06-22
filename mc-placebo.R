###################################
# Placebo MC estimates #
###################################

## Loading Source files
library(MCPanel)
library(glmnet)
library(ggplot2)

# Setup parallel processing 
library(parallel)
library(doParallel)

cores <- 4#detectCores()

cl <- parallel::makeForkCluster(cores)

doParallel::registerDoParallel(cores) # register cores (<p)

RNGkind("L'Ecuyer-CMRG") # ensure random number generation

outcome.vars <- c("CBWbord","CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_0","seekdur_1_2","seekdur_3more")

for(o in outcome.vars){
  print(o)
  
  ## Analysis 1: ST vs AT (retrospective, X=CBW) 
  
  # Eastern cluster
  print(paste0("Estimates for Analysis 1, Eastern cluster, outcome:",o))

  outcomes.cbw.eastern <- readRDS(paste0("data/outcomes-cbw-eastern-",o,".rds"))
  
  outcomes.cbw.eastern$mask[rownames(outcomes.cbw.eastern$mask)%in%outcomes.cbw.eastern$treated,] <- abs(outcomes.cbw.eastern$mask[rownames(outcomes.cbw.eastern$mask)%in%outcomes.cbw.eastern$treated,]-1) # estimate Y(1)_LT,pre
  
  # Discard post-treatment periods
  outcomes.cbw.eastern.placebo <- outcomes.cbw.eastern
  outcomes.cbw.eastern.placebo$M <- outcomes.cbw.eastern$M[,1:which(colnames(outcomes.cbw.eastern$M)=="20081")-1]
  outcomes.cbw.eastern.placebo$mask <- outcomes.cbw.eastern$mask[,1:which(colnames(outcomes.cbw.eastern$mask)=="20081")-1]
  
  # Get p-values
  source("MCEst.R")
  source("ChernoTest.R")
  
  t_final_placebo <- ncol(outcomes.cbw.eastern.placebo$M ) # all periods 
  
  taus <- 1:length((4:t_final_placebo))
  
  treat_indices_order <- outcomes.cbw.eastern.placebo$treated
  
  moving.block.placebo <- foreach(t = taus) %dopar% {
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.cbw.eastern.placebo[c("M","mask")], ns=1000, treat_indices_order=treat_indices_order, permtype="moving.block",t0=t0_placebo,rev=TRUE,covars=FALSE)}
  saveRDS(moving.block.placebo,paste0("results/no-covars/moving-block-placebo-cbw-eastern-",o,".rds"))
  
  iid.block.placebo <- foreach(t = taus) %dopar% {
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.cbw.eastern.placebo[c("M","mask")], ns=1000, treat_indices_order=treat_indices_order, permtype="iid.block",t0=t0_placebo,rev=TRUE,covars=FALSE)}
  saveRDS(iid.block.placebo,paste0("results/no-covars/iid-block-placebo-cbw-eastern-",o,".rds"))
  
  iid.placebo <- foreach(t = taus) %dopar% {
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.cbw.eastern.placebo[c("M","mask")],ns=1000, treat_indices_order=treat_indices_order, permtype="iid",t0=t0_placebo,rev=TRUE,covars=FALSE)}
  saveRDS(iid.placebo,paste0("results/no-covars/iid-placebo-cbw-eastern-",o,".rds"))
  
  # Swiss cluster
  print(paste0("Estimates for Analysis 1, swiss cluster, outcome:",o))
  
  outcomes.cbw.swiss <- readRDS(paste0("data/outcomes-cbw-swiss-",o,".rds"))
  
  outcomes.cbw.swiss$mask[rownames(outcomes.cbw.swiss$mask)%in%outcomes.cbw.swiss$treated,] <- abs(outcomes.cbw.swiss$mask[rownames(outcomes.cbw.swiss$mask)%in%outcomes.cbw.swiss$treated,]-1) # estimate Y(1)_LT,pre
  
  # Discard post-treatment periods
  outcomes.cbw.swiss.placebo <- outcomes.cbw.swiss
  outcomes.cbw.swiss.placebo$M <- outcomes.cbw.swiss$M[,1:which(colnames(outcomes.cbw.swiss$M)=="20072")-1]
  outcomes.cbw.swiss.placebo$mask <- outcomes.cbw.swiss$mask[,1:which(colnames(outcomes.cbw.swiss$mask)=="20072")-1]
  
  # Get p-values

  t_final_placebo <- ncol(outcomes.cbw.swiss.placebo$M ) # all periods 
  
  taus <- 1:length((4:t_final_placebo))
  
  treat_indices_order <- outcomes.cbw.swiss.placebo$treated
  
  moving.block.placebo <- foreach(t = taus) %dopar% {
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.cbw.swiss.placebo[c("M","mask")], ns=1000, treat_indices_order=treat_indices_order, permtype="moving.block",t0=t0_placebo,rev=TRUE,covars=FALSE)}
  saveRDS(moving.block.placebo,paste0("results/no-covars/moving-block-placebo-cbw-swiss-",o,".rds"))
  
  iid.block.placebo <- foreach(t = taus) %dopar% {
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.cbw.swiss.placebo[c("M","mask")], ns=1000, treat_indices_order=treat_indices_order, permtype="iid.block",t0=t0_placebo,rev=TRUE,covars=FALSE)}
  saveRDS(iid.block.placebo,paste0("results/no-covars/iid-block-placebo-cbw-swiss-",o,".rds"))
  
  iid.placebo <- foreach(t = taus) %dopar% {
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.cbw.swiss.placebo[c("M","mask")],ns=1000, treat_indices_order=treat_indices_order, permtype="iid",t0=t0_placebo,rev=TRUE,covars=FALSE)}
  saveRDS(iid.placebo,paste0("results/no-covars/iid-placebo-cbw-swiss-",o,".rds"))
  
  ## Analysis 2: ST vs NT (forward, X=LM)
  
  # Eastern cluster
  print(paste0("Estimates for Analysis 1, Eastern cluster, outcome:",o))
  
  outcomes.lm.eastern <- readRDS(paste0("data/outcomes-lm-eastern-",o,".rds"))
  
  # Discard post-treatment periods
  outcomes.lm.eastern.placebo <- outcomes.lm.eastern
  outcomes.lm.eastern.placebo$M <- outcomes.lm.eastern$M[,1:which(colnames(outcomes.lm.eastern$M)=="20081")-1]
  outcomes.lm.eastern.placebo$mask <- outcomes.lm.eastern$mask[,1:which(colnames(outcomes.lm.eastern$mask)=="20081")-1]
  
  # Get p-values

  t_final_placebo <- ncol(outcomes.lm.eastern.placebo$M ) # all periods 
  
  taus <- 1:length((4:t_final_placebo))
  
  treat_indices_order <- outcomes.lm.eastern.placebo$treated
  
  moving.block.placebo <- foreach(t = taus) %dopar% {
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.lm.eastern.placebo[c("M","mask")], ns=1000, treat_indices_order=treat_indices_order, permtype="moving.block",t0=t0_placebo,rev=TRUE,covars=FALSE)}
  saveRDS(moving.block.placebo,paste0("results/no-covars/moving-block-placebo-lm-eastern-",o,".rds"))
  
  iid.block.placebo <- foreach(t = taus) %dopar% {
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.lm.eastern.placebo[c("M","mask")], ns=1000, treat_indices_order=treat_indices_order, permtype="iid.block",t0=t0_placebo,rev=TRUE,covars=FALSE)}
  saveRDS(iid.block.placebo,paste0("results/no-covars/iid-block-placebo-lm-eastern-",o,".rds"))
  
  iid.placebo <- foreach(t = taus) %dopar% {
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.lm.eastern.placebo[c("M","mask")],ns=1000, treat_indices_order=treat_indices_order, permtype="iid",t0=t0_placebo,rev=TRUE,covars=FALSE)}
  saveRDS(iid.placebo,paste0("results/no-covars/iid-placebo-lm-eastern-",o,".rds"))
  
  # Swiss cluster
  print(paste0("Estimates for Analysis 1, swiss cluster, outcome:",o))
  
  outcomes.lm.swiss <- readRDS(paste0("data/outcomes-lm-swiss-",o,".rds"))
  
  # Discard post-treatment periods
  outcomes.lm.swiss.placebo <- outcomes.lm.swiss
  outcomes.lm.swiss.placebo$M <- outcomes.lm.swiss$M[,1:which(colnames(outcomes.lm.swiss$M)=="20072")-1]
  outcomes.lm.swiss.placebo$mask <- outcomes.lm.swiss$mask[,1:which(colnames(outcomes.lm.swiss$mask)=="20072")-1]
  
  # Get p-values
  
  t_final_placebo <- ncol(outcomes.lm.swiss.placebo$M ) # all periods 
  
  taus <- 1:length((4:t_final_placebo))
  
  treat_indices_order <- outcomes.lm.swiss.placebo$treated
  
  moving.block.placebo <- foreach(t = taus) %dopar% {
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.lm.swiss.placebo[c("M","mask")], ns=1000, treat_indices_order=treat_indices_order, permtype="moving.block",t0=t0_placebo,rev=TRUE,covars=FALSE)}
  saveRDS(moving.block.placebo,paste0("results/no-covars/moving-block-placebo-lm-swiss-",o,".rds"))
  
  iid.block.placebo <- foreach(t = taus) %dopar% {
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.lm.swiss.placebo[c("M","mask")], ns=1000, treat_indices_order=treat_indices_order, permtype="iid.block",t0=t0_placebo,rev=TRUE,covars=FALSE)}
  saveRDS(iid.block.placebo,paste0("results/no-covars/iid-block-placebo-lm-swiss-",o,".rds"))
  
  iid.placebo <- foreach(t = taus) %dopar% {
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    ChernoTest(outcomes=outcomes.lm.swiss.placebo[c("M","mask")],ns=1000, treat_indices_order=treat_indices_order, permtype="iid",t0=t0_placebo,rev=TRUE,covars=FALSE)}
  saveRDS(iid.placebo,paste0("results/no-covars/iid-placebo-lm-swiss-",o,".rds"))
}