###################################
# Placebo MC estimates #
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
  
  # Use pre-treatment for LT (no missing values under the null)
  outcomes.cbw.placebo <- outcomes.cbw
  outcomes.cbw.placebo$mask <- outcomes.cbw$mask[,1:(which(colnames(outcomes.cbw$mask)=="20091")-1)][!rownames(outcomes.cbw$mask) %in% outcomes.cbw$control,] -1   # all zeros 
  outcomes.cbw.placebo$M <- outcomes.cbw$M[,1:(which(colnames(outcomes.cbw$mask)=="20091")-1)][!rownames(outcomes.cbw$mask) %in% outcomes.cbw$control,]
  outcomes.cbw.placebo$W <- outcomes.cbw$W[,1:(which(colnames(outcomes.cbw$mask)=="20091")-1)][!rownames(outcomes.cbw$mask) %in% outcomes.cbw$control,]
  
  # Get optimal stationary bootstrap lengths
  source("PolitisWhite.R")
  
  bopt <- b.star(t(outcomes.cbw$M),round=TRUE)[,1]
  
  # Bootstrap for per-period effects
  source("MCEstBoot.R")
  
  boot <- tsboot(tseries=ts(t(outcomes.cbw.placebo$M)), MCEstBoot, mask=outcomes.cbw.placebo$mask, W=outcomes.cbw.placebo$W, covars=FALSE, rev=TRUE, R=999, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot, paste0("results/placebo-boot-cbw-",o,".rds")) 
}