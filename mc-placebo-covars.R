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
  
  # Use pre-treatment for LT (no missing values)
  outcomes.cbw.placebo <- outcomes.cbw
  outcomes.cbw.placebo$mask <- outcomes.cbw$mask[rownames(outcomes.cbw$mask)%in%outcomes.cbw$treated,][,1:(which(colnames(outcomes.cbw$mask)=="20091")-1)] -1 # all zeros
  outcomes.cbw.placebo$M <- outcomes.cbw$M[,colnames(outcomes.cbw$M)%in%colnames(outcomes.cbw.placebo$mask)][rownames(outcomes.cbw$M)%in%rownames(outcomes.cbw.placebo$mask),]
  outcomes.cbw.placebo$W <- outcomes.cbw$W[,colnames(outcomes.cbw$W)%in%colnames(outcomes.cbw.placebo$mask)][rownames(outcomes.cbw$W)%in%rownames(outcomes.cbw.placebo$mask),]
  outcomes.cbw.placebo$X <- outcomes.cbw$X[,colnames(outcomes.cbw$X)%in%colnames(outcomes.cbw.placebo$mask)][rownames(outcomes.cbw$X)%in%rownames(outcomes.cbw.placebo$mask),]
  outcomes.cbw.placebo$X.hat <- outcomes.cbw$X.hat[,colnames(outcomes.cbw$X.hat)%in%colnames(outcomes.cbw.placebo$mask)][rownames(outcomes.cbw$X.hat)%in%rownames(outcomes.cbw.placebo$mask),]
  
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
    t0_placebo <- t_final_placebo-t 
    tsboot(tseries=t(outcomes.cbw.placebo$M), MCEstBoot, mask=outcomes.cbw.placebo$mask, W=outcomes.cbw.placebo$W, X=outcomes.cbw.placebo$X,X.hat=outcomes.cbw.placebo$X.hat, eastern=outcomes.cbw$eastern, covars=TRUE, rev=TRUE, t0=t0_placebo, R=1000, parallel = "multicore", l=bopt, sim = "geom")})
  saveRDS(boot.trajectory.eastern.placebo.cbw,paste0("results/boot-trajectory-eastern-placebo-cbw-",o,"-covars.rds"))
  
  boot.trajectory.swiss.placebo.cbw <- lapply(taus, function(t){
    t0_placebo <- t_final_placebo-t # 
    tsboot(tseries=t(outcomes.cbw.placebo$M), MCEstBoot, mask=outcomes.cbw.placebo$mask, W=outcomes.cbw.placebo$W, X=outcomes.cbw.placebo$X,X.hat=outcomes.cbw.placebo$X.hat, swiss=outcomes.cbw$swiss, covars=TRUE, rev=TRUE, t0=t0_placebo, R=1000, parallel = "multicore", l=bopt, sim = "geom")})
  saveRDS(boot.trajectory.swiss.placebo.cbw,paste0("results/boot-trajectory-swiss-placebo-cbw-",o,"-covars.rds"))
  
  ## Analysis 2: ST vs NT (forward, X=LM)
  
  print(paste0("Estimates for Analysis 1, Eastern cluster, outcome:",o))
  
  outcomes.lm <- readRDS(paste0("data/outcomes-lm-",o,".rds"))
  
  # Use pre-treatment (no missing values)
  outcomes.lm.placebo <- outcomes.lm
  outcomes.lm.placebo$mask <- outcomes.lm$mask[,1:(which(colnames(outcomes.lm$mask)=="20072")-1)] # all zeros
  outcomes.lm.placebo$M <- outcomes.lm$M[,colnames(outcomes.lm$M)%in%colnames(outcomes.lm.placebo$mask)][rownames(outcomes.lm$M)%in%rownames(outcomes.lm.placebo$mask),]
  outcomes.lm.placebo$W <- outcomes.lm$W[,colnames(outcomes.lm$W)%in%colnames(outcomes.lm.placebo$mask)][rownames(outcomes.lm$W)%in%rownames(outcomes.lm.placebo$mask),]
  outcomes.lm.placebo$X <- outcomes.lm$X[,colnames(outcomes.lm$X)%in%colnames(outcomes.lm.placebo$mask)][rownames(outcomes.lm$X)%in%rownames(outcomes.lm.placebo$mask),]
  outcomes.lm.placebo$X.hat <- outcomes.lm$X.hat[,colnames(outcomes.lm$X.hat)%in%colnames(outcomes.lm.placebo$mask)][rownames(outcomes.lm$X.hat)%in%rownames(outcomes.lm.placebo$mask),]
  
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