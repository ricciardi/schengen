###################################
#  MC estimates #
###################################

## Libraries
library(MCPanel)
library(boot)

# Setup parallel processing 
library(parallel)
library(doParallel)

cores <- detectCores()

cl <- parallel::makeForkCluster(cores)

doParallel::registerDoParallel(cores) # register cores (<p)

RNGkind("L'Ecuyer-CMRG") # ensure random number generation

outcome.vars <- c("CBWbord","CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_3more")

for(o in outcome.vars){
  print(o)

  ## Analysis 1: ST vs AT (retrospective, X=CBW) 
  
  print(paste0("Estimates for Analysis 1, outcome:",o))
  
  outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  # Get treatment effect estimates
    
  source('MCEst.R')
  mc.estimates.cbw <- MCEst(outcomes.cbw, rev=TRUE, covars=FALSE)
  saveRDS(mc.estimates.cbw, paste0("results/mc-estimates-cbw-",o,".rds"))
  
  # Get optimal stationary bootstrap lengths
  source("PolitisWhite.R")
  
  bopt <- b.star(t(outcomes.cbw$M),round=TRUE)[,1]
  
  # Bootstrap for per-period effects
  source("MCEstBoot.R")
  
  t0.eastern <- which(colnames(outcomes.cbw$mask)==20111)
  t0.swiss <- which(colnames(outcomes.cbw$mask)==20072)
  
  boot <- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W,covars=FALSE, rev=TRUE, R=999, parallel = "multicore", l=bopt, sim = "geom") 
  saveRDS(boot, paste0("results/boot-cbw-",o,".rds")) 
  
  # boot.eastern <- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, t0=t0.eastern, eastern=outcomes.cbw$eastern, covars=FALSE, rev=TRUE, R=999, parallel = "multicore", l=bopt, sim = "geom") 
  # saveRDS(boot.eastern, paste0("results/boot-cbw-eastern",o,".rds")) 
  # 
  # boot.swiss <- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, t0=t0.swiss, swiss=outcomes.cbw$swiss, covars=FALSE, rev=TRUE, R=999, parallel = "multicore", l=bopt, sim = "geom") 
  # saveRDS(boot.swiss, paste0("results/boot-cbw-swiss",o,".rds")) 

  # Bootstrap for trajectories
  # Resample trajectories without time component, calculate ATTs for each cluster
  source("MCEstBootTraj.R")
  
  impact <- mc.estimates.cbw$impact # = boot_result$t0
  
  trajectory.eastern <- rowMeans(impact[,1:(t0.eastern-1)])
  trajectory.swiss <- rowMeans(impact[,1:(t0.swiss-1)])
  
  # eastern

  boot.trajectory.eastern <- boot(trajectory.eastern, 
                                  MCEstBootTraj, 
                                  eastern=outcomes.cbw$eastern,
                                  R=999,
                                  parallel = "multicore") 
  
  print(boot.trajectory.eastern$t0)
  print(boot.ci(boot.trajectory.eastern,type=c("norm","basic", "perc")))
  
  saveRDS(boot.trajectory.eastern, paste0("results/boot-cbw-trajectory-eastern-",o,".rds")) 
  
  # swiss
  
  boot.trajectory.swiss <- boot(trajectory.swiss, 
                                  MCEstBootTraj, 
                                  swiss=outcomes.cbw$swiss,
                                  R=999,
                                  parallel = "multicore") 
  
  print(boot.trajectory.swiss$t0)
  print(boot.ci(boot.trajectory.swiss,type=c("norm","basic", "perc")))
  
  saveRDS(boot.trajectory.swiss, paste0("results/boot-cbw-trajectory-swiss-",o,".rds")) 
}