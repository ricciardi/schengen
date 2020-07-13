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
  
  # Use post-treatment (all zeros)
  outcomes.cbw.placebo <- outcomes.cbw
  outcomes.cbw.placebo$mask <- outcomes.cbw$mask[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  outcomes.cbw.placebo$M <- outcomes.cbw$M[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  outcomes.cbw.placebo$W <- outcomes.cbw$W[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  
  # Get optimal stationary bootstrap lengths
  source("PolitisWhite.R")
  
  bopt <- b.star(t(outcomes.cbw$M),round=TRUE)[,1]
  
  # Random staggered adoption among actual treated 
  T0 <- round(c(ncol(outcomes.cbw.placebo$mask)-1, ncol(outcomes.cbw.placebo$mask)/1.25, ncol(outcomes.cbw.placebo$mask)/1.5)) # vary t0
  boot <- lapply(T0, function(t0){
    treat_indices <- which(rownames(outcomes.cbw.placebo$mask) %in%outcomes.cbw.placebo$treated) # keep treated fixed to actual treated
    treat_mat <- (1-stag_adapt(outcomes.cbw.placebo$M, length(treat_indices),t0, treat_indices)) # invert again in MCEst
    
    rotate <- function(x) t(apply(x, 2, rev))
    
    outcomes.cbw.placebo$mask <- rotate(rotate(treat_mat)) # retrospective analysis
    
    source('MCEst.R')
    mc.estimates.cbw.placebo <- MCEst(outcomes.cbw.placebo, rev=TRUE, covars=FALSE)
    
    # Resample trajectories without time component, calculate ATTs for each cluster
    source("MCEstBootTraj.R")
    
    impact <- mc.estimates.cbw.placebo$impact # = boot_result$t0
    
    trajectory.eastern <- rowMeans(impact[,1:(t0-1)])
    trajectory.swiss <- rowMeans(impact[,1:(t0-1)])
    
    # eastern
    
    boot.trajectory.eastern <- boot(trajectory.eastern, 
                                    MCEstBootTraj, 
                                    eastern=outcomes.cbw$eastern,
                                    R=999,
                                    parallel = "multicore") 
    
  
    # swiss
    
    boot.trajectory.swiss <- boot(trajectory.swiss, 
                                  MCEstBootTraj, 
                                  swiss=outcomes.cbw$swiss,
                                  R=999,
                                  parallel = "multicore") 
    
    return(list("eastern"=boot.trajectory.eastern,"swiss"=boot.trajectory.swiss))

  })
  names(boot) <- T0/ncol(outcomes.cbw.placebo$mask)
  saveRDS(boot, paste0("results/placebo-boot-cbw-",o,".rds")) 
}