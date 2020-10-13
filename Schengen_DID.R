###################################
#  DID estimates #
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

outcome.vars <- c("CBWbord","CBWbordEMPL","Thwusual")

for(o in outcome.vars){
  print(paste0("Outcome: ", o))

  print(paste0("No covariates + FEs, outcome:",o))
  
  outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  
  # Get optimal stationary bootstrap lengths
  source("PolitisWhite.R")
  
  bopt <- b.star(t(outcomes.cbw$M),round=TRUE)[,1]
  
  # Bootstrap for per-period effects
  source("DIDEstBoot.R")
  
  boot <- tsboot(tseries=ts(t(outcomes.cbw$M)), DIDEstBoot, mask=outcomes.cbw$mask, rev=TRUE, R=999, parallel = "multicore", l=bopt, sim = "geom") 
  
  # Bootstrap for trajectories
  # Resample trajectories without time component, calculate ATTs for each cluster
  source("MCEstBootTraj.R")
  
  impact <- boot$t0

  t0.eastern <- which(colnames(outcomes.cbw$mask)==20111)
  t0.swiss <- which(colnames(outcomes.cbw$mask)==20091)
  
  trajectory.eastern <- rowMeans(impact[,1:(t0.eastern-1)]) # Schengen + FoM
  trajectory.swiss <- rowMeans(impact[,1:(t0.swiss-1)]) # Schengen + FoM

  # eastern
  
  boot.trajectory.eastern <- boot(trajectory.eastern, 
                                  MCEstBootTraj, 
                                  eastern=outcomes.cbw$eastern,
                                  R=999,
                                  parallel = "multicore") 
  
  print(paste0("Eastern: Combined treatment effect: ", boot.trajectory.eastern$t0))
  print(boot.ci(boot.trajectory.eastern,type=c("norm","basic", "perc")))
  
  # swiss
  
  boot.trajectory.swiss <- boot(trajectory.swiss, 
                                MCEstBootTraj, 
                                swiss=outcomes.cbw$swiss,
                                R=999,
                                parallel = "multicore") 
  
  print(paste0("Swiss: Combined treatment effect: ", boot.trajectory.swiss$t0))
  print(boot.ci(boot.trajectory.swiss,type=c("norm","basic", "perc")))
}