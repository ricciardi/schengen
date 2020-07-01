###################################
# Placebo MC estimates with covariates #
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

outcome.vars <- c("N_CBWbord","CBWbord","CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_0","seekdur_1_2","seekdur_3more")

for(o in outcome.vars){
  
  ## Analysis 1: ST vs AT (retrospective, X=CBW) 
  
  print(paste0("Estimates for Analysis 1, outcome:",o))
  
  outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  # Discard pre-treatment periods (retrospective)
  outcomes.cbw.placebo <- outcomes.cbw
  outcomes.cbw.placebo$M <- outcomes.cbw$M[,which(colnames(outcomes.cbw$M)=="20091"):ncol(outcomes.cbw$M)] 
  outcomes.cbw.placebo$W <- outcomes.cbw$W[,which(colnames(outcomes.cbw$W)=="20091"):ncol(outcomes.cbw$M)]
  outcomes.cbw.placebo$mask <- outcomes.cbw$mask[,which(colnames(outcomes.cbw$mask)=="20091"):ncol(outcomes.cbw$M)]
  outcomes.cbw.placebo$X <- outcomes.cbw$X[,which(colnames(outcomes.cbw$X)=="20091"):ncol(outcomes.cbw$M)]
  outcomes.cbw.placebo$X.hat <- outcomes.cbw$X.hat[,which(colnames(outcomes.cbw$X.hat)=="20091"):ncol(outcomes.cbw$M)]
  
  # Get optimal stationary bootstrap lengths
  source("PolitisWhite.R")
  
  bopt <- try(b.star(t(outcomes.cbw.placebo$M),round=TRUE)[,1])  # get optimal bootstrap lengths
  
  if("try-error" %in% class(bopt)) bopt <- rep(3,ncol(outcomes.cbw.placebo$M))
  
  # Get p-values
  source("MCEst.R")
  source("MCEstBoot.R")
  
  t_final_placebo <- ncol(outcomes.cbw.placebo$M ) # all periods 
  
  taus <- 1:5
  
  boot.trajectory.eastern.placebo.cbw <- lapply(taus, function(t){
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    tsboot(tseries=t(outcomes.cbw.placebo$M), MCEstBoot, mask=outcomes.cbw.placebo$mask, W=outcomes.cbw.placebo$W, X=outcomes.cbw.placebo$X,X.hat=outcomes.cbw.placebo$X.hat, eastern=outcomes.cbw$eastern, covars=TRUE, rev=TRUE, t0=t0_placebo, R=1000, parallel = "multicore", l=bopt, sim = "geom")})
  saveRDS(boot.trajectory.eastern.placebo.cbw,paste0("results/boot-trajectory-eastern-placebo-cbw-",o,"-covars.rds"))
  
  boot.trajectory.swiss.placebo.cbw <- lapply(taus, function(t){
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    tsboot(tseries=t(outcomes.cbw.placebo$M), MCEstBoot, mask=outcomes.cbw.placebo$mask, W=outcomes.cbw.placebo$W, X=outcomes.cbw.placebo$X,X.hat=outcomes.cbw.placebo$X.hat, swiss=outcomes.cbw$swiss, covars=TRUE, rev=TRUE, t0=t0_placebo, R=1000, parallel = "multicore", l=bopt, sim = "geom")})
  saveRDS(boot.trajectory.swiss.placebo.cbw,paste0("results/boot-trajectory-swiss-placebo-cbw-",o,"-covars.rds"))
  
  ## Analysis 2: ST vs NT (forward, X=LM)
  
  print(paste0("Estimates for Analysis 1, Eastern cluster, outcome:",o))
  
  outcomes.lm <- readRDS(paste0("data/outcomes-lm-",o,".rds"))
  
  # Discard post-treatment periods
  outcomes.lm.placebo <- outcomes.lm
  outcomes.lm.placebo$M <- outcomes.lm$M[,1:which(colnames(outcomes.lm$M)=="20091")-1]
  outcomes.lm.placebo$W <- outcomes.lm$W[,1:which(colnames(outcomes.lm$W)=="20091")-1]
  outcomes.lm.placebo$mask <- outcomes.lm$mask[,1:which(colnames(outcomes.lm$mask)=="20091")-1]
  outcomes.lm.placebo$X <- outcomes.lm$X[,which(colnames(outcomes.lm$X)=="20091"):ncol(outcomes.lm$M)]
  outcomes.lm.placebo$X.hat <- outcomes.lm$X.hat[,which(colnames(outcomes.lm$X.hat)=="20091"):ncol(outcomes.lm$M)]
  
  # Get p-values
  
  t_final_placebo <- ncol(outcomes.lm.placebo$M ) # all periods 
  
  taus <- 1:5
  
  boot.trajectory.eastern.placebo.lm <- lapply(taus, function(t){
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    tsboot(tseries=t(outcomes.lm.placebo$M), MCEstBoot, mask=outcomes.lm.placebo$mask, W=outcomes.lm.placebo$W,X=outcomes.lm.placebo$X,X.hat=outcomes.lm.placebo$X.hat, eastern=outcomes.lm$eastern, covars=TRUE, rev=FALSE, t0=t0_placebo, R=1000, parallel = "multicore", l=bopt, sim = "geom")})
  saveRDS(boot.trajectory.eastern.placebo.lm,paste0("results/boot-trajectory-eastern-placebo-lm-",o,"-covars.rds"))
  
  boot.trajectory.swiss.placebo.lm <- lapply(taus, function(t){
    t0_placebo <- t_final_placebo-t # n pre-treatment periods
    tsboot(tseries=t(outcomes.lm.placebo$M), MCEstBoot, mask=outcomes.lm.placebo$mask, W=outcomes.lm.placebo$W,X=outcomes.lm.placebo$X,X.hat=outcomes.lm.placebo$X.hat, swiss=outcomes.lm$swiss, covars=TRUE, rev=FALSE, t0=t0_placebo, R=1000, parallel = "multicore", l=bopt, sim = "geom")})
  saveRDS(boot.trajectory.swiss.placebo.lm,paste0("results/boot-trajectory-swiss-placebo-lm-",o,"-covars.rds"))
}