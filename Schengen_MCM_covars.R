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
  
  # Bootstrap for trajectories
  # Resample trajectories without time component, calculate ATTs for each cluster
  source("MCEstBootTraj.R")
  
  impact.eastern <- mc.estimates.cbw.eastern$impact 
  impact.swiss <- mc.estimates.cbw.swiss$impact

  t0.eastern <- which(colnames(outcomes.cbw$mask)==20111) # t0-1 in MCEstBootTraj
  t0.swiss <- which(colnames(outcomes.cbw$mask)==20091)

  # eastern
  
  boot.trajectory.eastern <- boot(impact.eastern, 
                                  MCEstBootTraj, 
                                  R=999,
                                  t0.eastern=t0.eastern,
                                  eastern=outcomes.cbw$eastern,
                                  parallel = "multicore") 
  
  print(paste0("Eastern: Combined treatment effect (20051-20104): ", boot.trajectory.eastern$t0))
  print(boot.ci(boot.trajectory.eastern,type="perc"))
  
  cvs.eastern <- quantile(boot.trajectory.eastern$t0 - boot.trajectory.eastern$t + 1, c(0.025, 0.975))
  print(paste0("Eastern: Combined treatment effect (20051-20104): ", mean(boot.trajectory.eastern$t > cvs.eastern[1] & boot.trajectory.eastern$t < cvs.eastern[2]))) # p-val for percentile bootstrap
  
  saveRDS(boot.trajectory.eastern, paste0("results/boot-cbw-trajectory-eastern-",o,"-covars.rds")) 

  # swiss
  
  boot.trajectory.swiss <- boot(impact.swiss, 
                                MCEstBootTraj, 
                                t0.swiss=t0.swiss,
                                swiss=outcomes.cbw$swiss,
                                R=999,
                                parallel = "multicore") 
  
  print(paste0("Swiss: Combined treatment effect (20051-20084): ", boot.trajectory.swiss$t0))
  print(boot.ci(boot.trajectory.swiss,type="perc"))
  
  cvs.swiss <- quantile(boot.trajectory.swiss$t0 - boot.trajectory.swiss$t + 1, c(0.025, 0.975))
  print(paste0("swiss: Combined treatment effect (20051-20084): ", mean(boot.trajectory.swiss$t > cvs.swiss[1] & boot.trajectory.swiss$t < cvs.swiss[2]))) # p-val for percentile bootstrap
  
  saveRDS(boot.trajectory.swiss, paste0("results/boot-cbw-trajectory-swiss-",o,"-covars.rds")) 
  
  # eastern (1)
  
  boot.trajectory.eastern.1 <- boot(impact.eastern, 
                                  MCEstBootTraj, 
                                  R=999,
                                  t0.eastern=which(colnames(outcomes.cbw$mask)==20081),
                                  eastern=outcomes.cbw$eastern,
                                  parallel = "multicore") 
  
  print(paste0("Eastern: partial treatment effect (20051-20074): ", boot.trajectory.eastern.1$t0))
  print(boot.ci(boot.trajectory.eastern.1,type="perc"))
  
  cvs.eastern.1 <- quantile(boot.trajectory.eastern.1$t0 - boot.trajectory.eastern.1$t + 1, c(0.025, 0.975))
  print(paste0("Eastern: Combined treatment effect (20051-20074): ", mean(boot.trajectory.eastern.1$t > cvs.eastern.1[1] & boot.trajectory.eastern.1$t < cvs.eastern.1[2]))) # p-val for percentile bootstrap
  
  saveRDS(boot.trajectory.eastern.1, paste0("results/boot-cbw-trajectory-eastern-1-",o,"-covars.rds")) 
  
  # swiss (1)
  
  boot.trajectory.swiss.1 <- boot(impact.swiss, 
                                MCEstBootTraj, 
                                t0.swiss=which(colnames(outcomes.cbw$mask)==20071),
                                swiss=outcomes.cbw$swiss,
                                R=999,
                                parallel = "multicore") 
  
  print(paste0("Swiss: partial treatment effect (20051-20072): ", boot.trajectory.swiss.1$t0))
  print(boot.ci(boot.trajectory.swiss.1,type="perc"))
  
  cvs.swiss.1 <- quantile(boot.trajectory.swiss.1$t0 - boot.trajectory.swiss.1$t + 1, c(0.025, 0.975))
  print(paste0("Swiss: Combined treatment effect (20051-20072): ", mean(boot.trajectory.swiss.1$t > cvs.swiss.1[1] & boot.trajectory.swiss.1$t < cvs.swiss.1[2]))) # p-val for percentile bootstrap
  
  saveRDS(boot.trajectory.swiss.1, paste0("results/boot-cbw-trajectory-swiss-1-",o,"-covars.rds"))
  
  # eastern (2)
  
  boot.trajectory.eastern.2 <- boot(impact.eastern, 
                                    MCEstBootTraj, 
                                    R=999,
                                    t0.eastern=t0.eastern,
                                    eastern=outcomes.cbw$eastern,
                                    start= which(colnames(outcomes.cbw$mask)==20081),
                                    parallel = "multicore") 
  
  print(paste0("Eastern: partial treatment effect (20081-20104): ", boot.trajectory.eastern.2$t0))
  print(boot.ci(boot.trajectory.eastern.2,type="perc"))
  
  cvs.eastern.2 <- quantile(boot.trajectory.eastern.2$t0 - boot.trajectory.eastern.2$t + 1, c(0.025, 0.975))
  print(paste0("Eastern: Combined treatment effect (20081-20104): ", mean(boot.trajectory.eastern.2$t > cvs.eastern.2[1] & boot.trajectory.eastern.2$t < cvs.eastern.2[2]))) # p-val for percentile bootstrap
  
  saveRDS(boot.trajectory.eastern.2, paste0("results/boot-cbw-trajectory-eastern-2-",o,"-covars.rds")) 
  
  # swiss (2)
  
  boot.trajectory.swiss.2 <- boot(impact.swiss, 
                                  MCEstBootTraj, 
                                  t0.swiss=t0.swiss,
                                  swiss=outcomes.cbw$swiss,
                                  start= which(colnames(outcomes.cbw$mask)==20073),
                                  R=999,
                                  parallel = "multicore") 
  
  print(paste0("Swiss: partial treatment effect (20073-20084): ", boot.trajectory.swiss.2$t0))
  print(boot.ci(boot.trajectory.swiss.2,type="perc"))
  
  cvs.swiss.2 <- quantile(boot.trajectory.swiss.2$t0 - boot.trajectory.swiss.2$t + 1, c(0.025, 0.975))
  print(paste0("Swiss: Combined treatment effect (20073-20084): ", mean(boot.trajectory.swiss.2$t > cvs.swiss.2[1] & boot.trajectory.swiss.2$t < cvs.swiss.2[2]))) # p-val for percentile bootstrap
  
  saveRDS(boot.trajectory.swiss.2, paste0("results/boot-cbw-trajectory-swiss-2-",o,"-covars.rds")) 
}