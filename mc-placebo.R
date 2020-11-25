###################################
# Placebo MC estimates #
###################################

## Loading Source files
library(MCPanel)
library(boot)
library(Matrix)

# Setup parallel processing
library(parallel)
library(doParallel)

cores <- detectCores()

cl <- parallel::makeForkCluster(cores)

doParallel::registerDoParallel(cores) # register cores (<p)

RNGkind("L'Ecuyer-CMRG") # ensure random number generation

outcome.vars <- c("CBWbord","CBWbordEMPL")
sim.labels <- c("Staggered adoption","Simultaneous adoption")

for(i in c(0,1)){
  for(o in outcome.vars){
    
    ## Analysis 1: ST vs AT (retrospective, X=CBW) 
    
    print(paste0("Estimates for Analysis 1, outcome:",o,sim.labels[i+1]))
    
    outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
    
    # Use post-treatment (all zeros)
    outcomes.cbw.placebo <- outcomes.cbw
    outcomes.cbw.placebo$mask <- outcomes.cbw$mask[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
    outcomes.cbw.placebo$M <- outcomes.cbw$M[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
    outcomes.cbw.placebo$W <- outcomes.cbw$W[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
    outcomes.cbw.placebo$X <- outcomes.cbw$X[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
    outcomes.cbw.placebo$X.hat <- outcomes.cbw$X.hat[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
    
    # Random staggered adoption among actual treated 
    T0 <- round(c(ncol(outcomes.cbw.placebo$mask)-1, ncol(outcomes.cbw.placebo$mask)/1.25, ncol(outcomes.cbw.placebo$mask)/1.5, ncol(outcomes.cbw.placebo$mask)/2))
    boot <- lapply(T0, function(t0){
      treat_indices <- which(rownames(outcomes.cbw.placebo$M) %in%outcomes.cbw.placebo$treated) # keep treated fixed to actual treated
      if(i==1){
        treat_mat <- (1-stag_adapt(outcomes.cbw.placebo$M, length(treat_indices),t0, treat_indices)) # invert again in MCEst
      } else{
        treat_mat <- (1-simul_adapt(outcomes.cbw.placebo$M, length(treat_indices),t0, treat_indices)) # invert again in MCEst
      }
      
      rotate <- function(x) t(apply(x, 2, rev))
      
      outcomes.cbw.placebo$mask <- rotate(rotate(treat_mat)) # retrospective analysis
      rownames(outcomes.cbw.placebo$mask) <- rownames(outcomes.cbw.placebo$M )
      
      z.cbw.eastern <- c(rev(SSlogis(1:t0, Asym = 1, xmid = 0.85, scal = 8)),
                         SSlogis(1:(ncol(outcomes.cbw.placebo$mask)-t0), Asym = 1, xmid = 0.85, scal = 8))
      
      outcomes.cbw.placebo$z.cbw.eastern <- z.cbw.eastern
      outcomes.cbw.placebo$z.cbw.swiss <- z.cbw.eastern
      
      source('MCEst.R')
      mc.estimates.cbw.eastern.placebo <- MCEst(outcomes.cbw.placebo, cluster='eastern', rev=TRUE, covars=FALSE) 
      mc.estimates.cbw.swiss.placebo <- MCEst(outcomes.cbw.placebo, cluster='swiss', rev=TRUE, covars=FALSE)
      
      # Resample trajectories without time component, calculate ATTs for each cluster
      source("MCEstBootTraj.R")
      
      impact.eastern <- mc.estimates.cbw.eastern.placebo$impact
      impact.swiss <- mc.estimates.cbw.swiss.placebo$impact
      
      # eastern
      
      boot.trajectory.eastern <- boot(impact.eastern, 
                                      MCEstBootTraj, 
                                      t0.eastern=t0,
                                      eastern=outcomes.cbw.placebo$eastern,
                                      R=999,
                                      parallel = "multicore") 
      
      print(boot.trajectory.eastern$t0)
      print(boot.ci(boot.trajectory.eastern,type=c("norm","basic", "perc")))
      
      # swiss
      
      boot.trajectory.swiss <- boot(impact.swiss, 
                                    MCEstBootTraj, 
                                    t0.swiss=t0,
                                    swiss=outcomes.cbw.placebo$swiss,
                                    R=999,
                                    parallel = "multicore") 
      
      print(boot.trajectory.swiss$t0)
      print(boot.ci(boot.trajectory.swiss,type=c("norm","basic", "perc")))
      
      return(list("eastern"=boot.trajectory.eastern,"swiss"=boot.trajectory.swiss))
      
    })
    names(boot) <- T0/ncol(outcomes.cbw.placebo$mask)
    saveRDS(boot, paste0("results/placebo-boot-cbw-",o,i,".rds")) 
  }
}