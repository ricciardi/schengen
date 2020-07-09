###################################
# Placebo MC estimates with covariates #
###################################

## Loading Source files
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
  
  ## Analysis 1: ST vs AT (retrospective, X=CBW) 
  
  print(paste0("Estimates for Analysis 1, outcome:",o))
  
  outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  # Use pre-treatment
  outcomes.cbw.placebo <- outcomes.cbw
  outcomes.cbw.placebo$mask <- outcomes.cbw$mask[,1:(which(colnames(outcomes.cbw$mask)=="20072")-1)]
  outcomes.cbw.placebo$mask[outcomes.cbw.placebo$mask>0] <- 0
  outcomes.cbw.placebo$mask[rownames(outcomes.cbw.placebo$mask)%in%outcomes.cbw.placebo$eastern,][,1:(which(colnames(outcomes.cbw.placebo$mask)=="20052"))] <- 1
  outcomes.cbw.placebo$mask[rownames(outcomes.cbw.placebo$mask)%in%outcomes.cbw.placebo$swiss,][,1:(which(colnames(outcomes.cbw.placebo$mask)=="20054"))] <- 1
  outcomes.cbw.placebo$M <- outcomes.cbw$M[,1:(which(colnames(outcomes.cbw$mask)=="20072")-1)]
  outcomes.cbw.placebo$W <- outcomes.cbw$W[,1:(which(colnames(outcomes.cbw$mask)=="20072")-1)]
  outcomes.cbw.placebo$X <- outcomes.cbw$X[,1:(which(colnames(outcomes.cbw$mask)=="20072")-1)]
  outcomes.cbw.placebo$X.hat <- outcomes.cbw$X.hat[,1:(which(colnames(outcomes.cbw$mask)=="20072")-1)]
  
  # Get optimal stationary bootstrap lengths
  source("PolitisWhite.R")
  
  bopt <- b.star(t(outcomes.cbw$M),round=TRUE)[,1]
  
  # Bootstrap for per-period effects
  source("MCEstBoot.R")
  
  boot <- tsboot(tseries=ts(t(outcomes.cbw.placebo$M)), MCEstBoot, mask=outcomes.cbw.placebo$mask, W=outcomes.cbw.placebo$W, X=outcomes.cbw.placebo$X,X.hat=outcomes.cbw.placebo$X.hat,covars=TRUE, rev=TRUE, R=999, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot, paste0("results/placebo-boot-cbw-",o,"-covars.rds")) 
}