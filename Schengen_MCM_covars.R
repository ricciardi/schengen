###################################
#  MC estimates with covariates#
###################################

## Libraries
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

for(o in outcome.vars){
  print(paste0("Outcome: ", o))
  
  print(paste0("Covariates + FEs, outcome:",o))
  
  outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  # Get treatment effect estimates
    
  source('MCEst.R')
  mc.estimates.cbw <- MCEst(outcomes.cbw, rev=TRUE, covars=TRUE)
  saveRDS(mc.estimates.cbw, paste0("results/mc-estimates-cbw-",o,"-covars.rds"))
  
  print(paste0("Rank of L: ", mc.estimates.cbw$rankL))

  # Get optimal stationary bootstrap lengths
  source("PolitisWhite.R")
  
  bopt <- b.star(t(outcomes.cbw$M),round=TRUE)[,1]
  
  # Bootstrap for per-period effects
  source("MCEstBoot.R")
  
  boot <- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X, X.hat=outcomes.cbw$X.hat, 
                 z.cbw.eastern =outcomes.cbw$z.cbw.eastern, z.cbw.swiss = outcomes.cbw$z.cbw.swiss, eastern=outcomes.cbw$eastern, swiss=outcomes.cbw$swiss, covars=TRUE, rev=TRUE, R=999, parallel = "multicore", l=bopt, sim = "geom", best_L=mc.estimates.cbw$best_lambda_L, best_B=mc.estimates.cbw$best_lambda_B) 
  saveRDS(boot, paste0("results/boot-cbw-",o,"-covars.rds")) 
  
  # Estimates with no propensity score weighting
  outcomes.cbw$W <- outcomes.cbw$W.equal # equal weighting
  
  mc.estimates.cbw.equal <- MCEst(outcomes.cbw, rev=TRUE, covars=TRUE)
  saveRDS(mc.estimates.cbw.equal, paste0("results/mc-estimates-cbw-equal-",o,"-covars.rds"))
  
  print(paste0("Rank of L: ", mc.estimates.cbw.equal$rankL))
  
  boot.equal <- tsboot(tseries=ts(t(outcomes.cbw$M)), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X, X.hat=outcomes.cbw$X.hat, 
                 z.cbw.eastern =outcomes.cbw$z.cbw.eastern, z.cbw.swiss = outcomes.cbw$z.cbw.swiss, eastern=outcomes.cbw$eastern, swiss=outcomes.cbw$swiss, covars=TRUE, rev=TRUE, R=999, parallel = "multicore", l=bopt, sim = "geom", best_L=mc.estimates.cbw.equal$best_lambda_L, best_B=mc.estimates.cbw.equal$best_lambda_B) 
  saveRDS(boot.equal, paste0("results/boot-cbw-equal-",o,"-covars.rds")) 
  
  # Bootstrap for trajectories
  # Resample trajectories without time component, calculate ATTs for each cluster
  source("MCEstBootTraj.R")
  
  impact <- mc.estimates.cbw$impact # = boot_result$t0

  t0.eastern <- which(colnames(outcomes.cbw$mask)==20111)
  t0.swiss <- which(colnames(outcomes.cbw$mask)==20091)

  # eastern
  
  boot.trajectory.eastern <- boot(impact, 
                                  MCEstBootTraj, 
                                  R=999,
                                  t0.eastern=t0.eastern,
                                  eastern=outcomes.cbw$eastern,
                                  parallel = "multicore") 
  
  print(paste0("Eastern: Combined treatment effect: ", boot.trajectory.eastern$t0))
  print(boot.ci(boot.trajectory.eastern,type=c("norm","basic", "perc")))
  
  saveRDS(boot.trajectory.eastern, paste0("results/boot-cbw-trajectory-eastern-",o,"-covars.rds")) 

  # swiss
  
  boot.trajectory.swiss <- boot(impact, 
                                MCEstBootTraj, 
                                t0.swiss=t0.swiss,
                                swiss=outcomes.cbw$swiss,
                                R=999,
                                parallel = "multicore") 
  
  print(paste0("Swiss: Combined treatment effect: ", boot.trajectory.swiss$t0))
  print(boot.ci(boot.trajectory.swiss,type=c("norm","basic", "perc")))
  
  saveRDS(boot.trajectory.swiss, paste0("results/boot-cbw-trajectory-swiss-",o,"-covars.rds")) 
  
  # Estimates without propensity weighting
  impact.equal <- mc.estimates.cbw.equal$impact
  
  # eastern
  
  boot.trajectory.eastern.equal <- boot(impact.equal, 
                                        MCEstBootTraj, 
                                        t0.eastern=t0.eastern,
                                        eastern=outcomes.cbw$eastern,
                                        R=999,
                                        parallel = "multicore") 
  
  print(paste0("Eastern: Combined treatment effect: ", boot.trajectory.eastern.equal$t0))
  print(boot.ci(boot.trajectory.eastern.equal,type=c("norm","basic", "perc")))
  
  saveRDS(boot.trajectory.eastern.equal, paste0("results/boot-cbw-trajectory-eastern-equal-",o,"-covars.rds")) 
  
  # swiss
  
  boot.trajectory.swiss.equal <- boot(impact.equal, 
                                      MCEstBootTraj, 
                                      t0.swiss=t0.swiss,
                                      swiss=outcomes.cbw$swiss, 
                                      R=999,
                                      parallel = "multicore") 
  
  print(paste0("Swiss: Combined treatment effect: ", boot.trajectory.swiss.equal$t0))
  print(boot.ci(boot.trajectory.swiss.equal,type=c("norm","basic", "perc")))
  
  saveRDS(boot.trajectory.swiss.equal, paste0("results/boot-cbw-trajectory-swiss-equal-",o,"-covars.rds")) 
}