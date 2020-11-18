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

outcome.vars <- c("CBWbord","CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_0","seekdur_1_2","seekdur_3more")

for(o in outcome.vars){
  print(paste0("Outcome: ", o))
  
  print(paste0("Covariates + FEs, outcome:",o))
  
  outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  # Get treatment effect estimates
    
  source('MCEst.R')
  mc.estimates.cbw.eastern <- MCEst(outcomes.cbw, cluster='eastern', rev=TRUE, covars=TRUE)
  saveRDS(mc.estimates.cbw.eastern, paste0("results/mc-estimates-cbw-eastern-",o,"-covars.rds"))
  
  mc.estimates.cbw.swiss <- MCEst(outcomes.cbw, cluster='swiss', rev=TRUE, covars=TRUE)
  saveRDS(mc.estimates.cbw.swiss, paste0("results/mc-estimates-cbw-swiss-",o,"-covars.rds"))
  
  print(paste0("Rank of L (Eastern): ", mc.estimates.cbw.eastern$rankL))
  print(paste0("Rank of L (Swiss): ", mc.estimates.cbw.swiss$rankL))

  # Get optimal stationary bootstrap lengths
  source("PolitisWhite.R")
  
  bopt.eastern <- b.star(t(outcomes.cbw$M[!rownames(outcomes.cbw$M)%in%outcomes.cbw$swiss,]),round=TRUE)[,1]
  bopt.swiss <- b.star(t(outcomes.cbw$M[!rownames(outcomes.cbw$M)%in%outcomes.cbw$eastern,]),round=TRUE)[,1]
  
  # Bootstrap for per-period effects
  source("MCEstBoot.R")
  
  boot.eastern <- tsboot(tseries=t(outcomes.cbw$M), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X, X.hat=outcomes.cbw$X.hat, 
                 z.cbw.eastern =outcomes.cbw$z.cbw.eastern, z.cbw.swiss = outcomes.cbw$z.cbw.swiss, eastern=outcomes.cbw$eastern, swiss=outcomes.cbw$swiss, est_eastern=TRUE, covars=TRUE, rev=TRUE, R=999, parallel = "multicore", l=bopt.eastern, sim = "geom", best_L=mc.estimates.cbw.eastern$best_lambda_L, best_B=mc.estimates.cbw.eastern$best_lambda_B) 
  saveRDS(boot.eastern, paste0("results/boot-cbw-eastern-",o,"-covars.rds")) 
  
  boot.swiss <- tsboot(tseries=t(outcomes.cbw$M), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X, X.hat=outcomes.cbw$X.hat, 
                         z.cbw.eastern =outcomes.cbw$z.cbw.eastern, z.cbw.swiss = outcomes.cbw$z.cbw.swiss, eastern=outcomes.cbw$eastern, swiss=outcomes.cbw$swiss, est_swiss=TRUE, covars=TRUE, rev=TRUE, R=999, parallel = "multicore", l=bopt.swiss, sim = "geom", best_L=mc.estimates.cbw.swiss$best_lambda_L, best_B=mc.estimates.cbw.swiss$best_lambda_B) 
  saveRDS(boot.swiss, paste0("results/boot-cbw-swiss-",o,"-covars.rds")) 
  
  # Estimates with no propensity score weighting
  outcomes.cbw$W <- outcomes.cbw$W.equal # equal weighting
  
  mc.estimates.cbw.eastern.equal <- MCEst(outcomes.cbw, cluster='eastern', rev=TRUE, covars=TRUE)
  saveRDS(mc.estimates.cbw.eastern.equal, paste0("results/mc-estimates-cbw-eastern-equal-",o,"-covars.rds"))
  
  print(paste0("Rank of L (Eastern) (No weighting): ", mc.estimates.cbw.eastern.equal$rankL))
  
  mc.estimates.cbw.swiss.equal <- MCEst(outcomes.cbw, cluster='swiss', rev=TRUE, covars=TRUE)
  saveRDS(mc.estimates.cbw.swiss.equal, paste0("results/mc-estimates-cbw-swiss-equal-",o,"-covars.rds"))
  
  print(paste0("Rank of L (Swiss) (No weighting): ", mc.estimates.cbw.swiss.equal$rankL))
  
  boot.eastern.equal <- tsboot(tseries=t(outcomes.cbw$M), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X, X.hat=outcomes.cbw$X.hat, 
                 z.cbw.eastern =outcomes.cbw$z.cbw.eastern, z.cbw.swiss = outcomes.cbw$z.cbw.swiss, eastern=outcomes.cbw$eastern, swiss=outcomes.cbw$swiss, est_eastern=TRUE, covars=TRUE, rev=TRUE, R=999, parallel = "multicore", l=bopt.eastern, sim = "geom", best_L=mc.estimates.cbw.eastern.equal$best_lambda_L, best_B=mc.estimates.cbw.eastern.equal$best_lambda_B) 
  saveRDS(boot.eastern.equal, paste0("results/boot-cbw-eastern-equal-",o,"-covars.rds")) 
  
  boot.swiss.equal <- tsboot(tseries=t(outcomes.cbw$M), MCEstBoot, mask=outcomes.cbw$mask, W=outcomes.cbw$W, X=outcomes.cbw$X, X.hat=outcomes.cbw$X.hat, 
                               z.cbw.eastern =outcomes.cbw$z.cbw.eastern, z.cbw.swiss = outcomes.cbw$z.cbw.swiss, eastern=outcomes.cbw$eastern, swiss=outcomes.cbw$swiss, est_swiss=TRUE, covars=TRUE, rev=TRUE, R=999, parallel = "multicore", l=bopt.swiss, sim = "geom", best_L=mc.estimates.cbw.swiss.equal$best_lambda_L, best_B=mc.estimates.cbw.swiss.equal$best_lambda_B) 
  saveRDS(boot.eastern.equal, paste0("results/boot-cbw-eastern-equal-",o,"-covars.rds")) 
  
  # Bootstrap for trajectories
  # Resample trajectories without time component, calculate ATTs for each cluster
  source("MCEstBootTraj.R")
  
  impact.eastern <- mc.estimates.cbw.eastern$impact 
  impact.swiss <- mc.estimates.cbw.swiss$impact

  t0.eastern <- which(colnames(outcomes.cbw$mask)==20111)
  t0.swiss <- which(colnames(outcomes.cbw$mask)==20091)

  # eastern
  
  boot.trajectory.eastern <- boot(impact.eastern, 
                                  MCEstBootTraj, 
                                  R=999,
                                  t0.eastern=t0.eastern,
                                  eastern=outcomes.cbw$eastern,
                                  parallel = "multicore") 
  
  print(paste0("Eastern: Combined treatment effect: ", boot.trajectory.eastern$t0))
  print(boot.ci(boot.trajectory.eastern,type=c("norm","basic", "perc")))
  
  saveRDS(boot.trajectory.eastern, paste0("results/boot-cbw-trajectory-eastern-",o,"-covars.rds")) 

  # swiss
  
  boot.trajectory.swiss <- boot(impact.swiss, 
                                MCEstBootTraj, 
                                t0.swiss=t0.swiss,
                                swiss=outcomes.cbw$swiss,
                                R=999,
                                parallel = "multicore") 
  
  print(paste0("Swiss: Combined treatment effect: ", boot.trajectory.swiss$t0))
  print(boot.ci(boot.trajectory.swiss,type=c("norm","basic", "perc")))
  
  saveRDS(boot.trajectory.swiss, paste0("results/boot-cbw-trajectory-swiss-",o,"-covars.rds")) 
  
  # Estimates without propensity weighting
  impact.eastern.equal <- mc.estimates.cbw.eastern.equal$impact
  impact.swiss.equal <- mc.estimates.cbw.swiss.equal$impact
  
  # eastern
  
  boot.trajectory.eastern.equal <- boot(impact.eastern.equal, 
                                        MCEstBootTraj, 
                                        t0.eastern=t0.eastern,
                                        eastern=outcomes.cbw$eastern,
                                        R=999,
                                        parallel = "multicore") 
  
  print(paste0("Eastern: Combined treatment effect: ", boot.trajectory.eastern.equal$t0))
  print(boot.ci(boot.trajectory.eastern.equal,type=c("norm","basic", "perc")))
  
  saveRDS(boot.trajectory.eastern.equal, paste0("results/boot-cbw-trajectory-eastern-equal-",o,"-covars.rds")) 
  
  # swiss
  
  boot.trajectory.swiss.equal <- boot(impact.swiss.equal, 
                                      MCEstBootTraj, 
                                      t0.swiss=t0.swiss,
                                      swiss=outcomes.cbw$swiss, 
                                      R=999,
                                      parallel = "multicore") 
  
  print(paste0("Swiss: Combined treatment effect: ", boot.trajectory.swiss.equal$t0))
  print(boot.ci(boot.trajectory.swiss.equal,type=c("norm","basic", "perc")))
  
  saveRDS(boot.trajectory.swiss.equal, paste0("results/boot-cbw-trajectory-swiss-equal-",o,"-covars.rds")) 
}