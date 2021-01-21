###################################
#  DID and SCM estimates #
###################################

## Libraries
library(MCPanel)
library(boot)
library(Matrix)
library(NNLM)

# Setup parallel processing 
library(parallel)
library(doParallel)

cores <- detectCores()

cl <- parallel::makeForkCluster(cores)

doParallel::registerDoParallel(cores) # register cores (<p)

RNGkind("L'Ecuyer-CMRG") # ensure random number generation

outcome.vars <- c("CBWbord","CBWbordEMPL")
estimators <- c("DID","ADH")

for(estimator in estimators){
  for(o in outcome.vars){
    print(paste0("Estimator: ", estimator))
    print(paste0("Outcome: ", o))
    
    outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
    
    # Get optimal stationary bootstrap lengths
    source("PolitisWhite.R")
    
    bopt.eastern <- b.star(t(outcomes.cbw$M[!rownames(outcomes.cbw$M)%in%outcomes.cbw$swiss,]),round=TRUE)[,1]
    bopt.swiss <- b.star(t(outcomes.cbw$M[!rownames(outcomes.cbw$M)%in%outcomes.cbw$eastern,]),round=TRUE)[,1]
    
    # Bootstrap for per-period effects
    source("DIDEstBoot.R")
    
    boot.eastern <- tsboot(tseries=t(outcomes.cbw$M), DIDEstBoot, mask=outcomes.cbw$mask, eastern=outcomes.cbw$eastern, swiss=outcomes.cbw$swiss, est_eastern=TRUE, rev=TRUE, estimator=estimator, R=999, parallel = "multicore", l=bopt.eastern, sim = "geom") 
    boot.swiss <- tsboot(tseries=t(outcomes.cbw$M), DIDEstBoot, mask=outcomes.cbw$mask, eastern=outcomes.cbw$eastern, swiss=outcomes.cbw$swiss, est_swiss=TRUE, rev=TRUE, estimator=estimator, R=999, parallel = "multicore", l=bopt.swiss, sim = "geom") 
    
    # Bootstrap for trajectories
    # Resample trajectories without time component, calculate ATTs for each cluster
    source("MCEstBootTraj.R")
    
    impact.eastern <- boot.eastern$t0
    impact.swiss <- boot.swiss$t0
    
    t0.eastern <- which(colnames(outcomes.cbw$mask)==20111)
    t0.swiss <- which(colnames(outcomes.cbw$mask)==20091)
    
    # eastern
    
    boot.trajectory.eastern <- boot(impact.eastern, 
                                    MCEstBootTraj, 
                                    t0.eastern=t0.eastern,
                                    eastern=outcomes.cbw$eastern,
                                    R=999,
                                    parallel = "multicore") 
    
    print(paste0("Eastern: Combined treatment effect: ", boot.trajectory.eastern$t0))
    print(boot.ci(boot.trajectory.eastern,type=c("norm","basic", "perc")))
    
    boot.t.eastern.null <- boot.trajectory.eastern$t - mean(boot.trajectory.eastern$t,na.rm = TRUE) # center around zero
    
    boot.trajectory.eastern.pval <- (1+sum( abs(boot.t.eastern.null) > abs(boot.trajectory.eastern$t0)))/(999+1)
    
    print(paste0("Eastern p-val: Combined treatment effect (20051-20104): ", boot.trajectory.eastern.pval))
    
    # swiss
    
    boot.trajectory.swiss <- boot(impact.swiss, 
                                  MCEstBootTraj, 
                                  t0.swiss=t0.swiss,
                                  swiss=outcomes.cbw$swiss,
                                  R=999,
                                  parallel = "multicore") 
    
    print(paste0("Swiss: Combined treatment effect: ", boot.trajectory.swiss$t0))
    print(boot.ci(boot.trajectory.swiss,type=c("norm","basic", "perc")))
    
    boot.t.swiss.null <- boot.trajectory.swiss$t - mean(boot.trajectory.swiss$t,na.rm = TRUE) # center around zero
    
    boot.trajectory.swiss.pval <- (1+sum( abs(boot.t.swiss.null) > abs(boot.trajectory.swiss$t0)))/(999+1)
    
    print(paste0("Swiss p-val: Combined treatment effect (20051-20084): ", boot.trajectory.swiss.pval)) # p-val for percentile bootstrap
  }
}