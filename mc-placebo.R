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
  
  # Use pre-treatment
  outcomes.cbw.placebo <- outcomes.cbw
  outcomes.cbw.placebo$mask <- outcomes.cbw$mask[,1:(which(colnames(outcomes.cbw$mask)=="20072")-1)]
  outcomes.cbw.placebo$mask[outcomes.cbw.placebo$mask>0] <- 0
  outcomes.cbw.placebo$M <- outcomes.cbw$M[,1:(which(colnames(outcomes.cbw$mask)=="20072")-1)]
  outcomes.cbw.placebo$W <- outcomes.cbw$W[,1:(which(colnames(outcomes.cbw$mask)=="20072")-1)]
  
  # Get optimal stationary bootstrap lengths
  source("PolitisWhite.R")
  
  bopt <- b.star(t(outcomes.cbw$M),round=TRUE)[,1]
  
  # Random staggered adoption among actual treated 
  T0 <- 4:(ncol(outcomes.cbw.placebo$mask)-1) # vary t0
  boot <- lapply(T0, function(t0){
    treat_indices <- which(rownames(outcomes.cbw.placebo$mask) %in%outcomes.cbw.placebo$treated) # keep treated fixed to actual treated
    treat_mat <- (1-stag_adapt(outcomes.cbw.placebo$M, length(treat_indices),t0, treat_indices)) # invert again in MCEstBoot
    
    rotate <- function(x) t(apply(x, 2, rev))
    
    treat_mat <- rotate(rotate(treat_mat)) # retrospective analysis
    
    # Bootstrap for per-period effects
    source("MCEstBoot.R")
    
    tsboot(tseries=ts(t(outcomes.cbw.placebo$M)), MCEstBoot, mask=outcomes.cbw.placebo$mask, W=outcomes.cbw.placebo$W, covars=FALSE, rev=TRUE, R=1999, parallel = "multicore", l=bopt, sim = "geom") 
  })
  names(boot) <- T0
  saveRDS(boot, paste0("results/placebo-boot-cbw-",o,".rds")) 
}