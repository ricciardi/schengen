###################################
# Placebo MC estimates #
###################################

## Libraries
library(MCPanel)
library(Matrix)
library(data.table)
library(reshape)
library(reshape2)

source('MCEst.R')
source('utils.R')

## Set random seed
set.seed(10, "L'Ecuyer-CMRG")

outcome.vars <- c("CBWbord","CBWbordEMPL")
sim.labels <- c("Staggered adoption")

for(o in outcome.vars){
  
  ## Analysis 1: ST vs AT (retrospective, X=CBW) 
  
  print(paste0("Estimates for Analysis 1, outcome:",o,sim.labels[1]))
  
  outcomes <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  # Use post-treatment (all zeros)
  outcomes.placebo <- outcomes
  outcomes.placebo$mask <- outcomes$mask[,which(colnames(outcomes$mask)=="20111"):ncol(outcomes$mask)]
  outcomes.placebo$M <- outcomes$M[,which(colnames(outcomes$mask)=="20111"):ncol(outcomes$mask)]
  outcomes.placebo$W <- outcomes$W[,which(colnames(outcomes$mask)=="20111"):ncol(outcomes$mask)]
  outcomes.placebo$z_weights <- outcomes$z_weights[,which(colnames(outcomes$mask)=="20111"):ncol(outcomes$mask)]
  outcomes.placebo$X <- outcomes$X[,which(colnames(outcomes$mask)=="20111"):ncol(outcomes$mask)]
  outcomes.placebo$X.hat <- outcomes$X.hat[,which(colnames(outcomes$mask)=="20111"):ncol(outcomes$mask)]
  
  # Random staggered adoption among actual treated 
  T0 <- ceiling(c(ncol(outcomes.placebo$mask)*0.5, ncol(outcomes.placebo$mask)*0.65, ncol(outcomes.placebo$mask)*0.8, ncol(outcomes.placebo$mask)*0.95))
  for(t0 in T0){
    treat_indices <- which(rownames(outcomes.placebo$M) %in%outcomes.placebo$treated) # keep treated fixed to actual treated
    
    treat_mat <- (1-stag_adapt(outcomes.placebo$M, length(treat_indices),t0, treat_indices)) # invert again in MCEst
    
    outcomes.placebo$mask <- treat_mat[,c(ncol(outcomes.placebo$M):1)]  # retrospective analysis
    rownames(outcomes.placebo$mask) <- rownames(outcomes.placebo$M )
    
    mc.estimates.eastern.placebo <- MCEst(outcomes.placebo, cluster='eastern', rev=TRUE, covars=TRUE)
    mc.estimates.swiss.placebo <- MCEst(outcomes.placebo, cluster='swiss', rev=TRUE, covars=TRUE)
    
    # eastern
    
    att.eastern <- apply(mc.estimates.eastern.placebo$impact[which(rownames(mc.estimates.eastern.placebo$impact) %in% outcomes.placebo$eastern),]*outcomes.placebo$mask[which(rownames(outcomes.placebo$mask) %in% outcomes.placebo$eastern),],1,nzmean)[outcomes.placebo$eastern] 
    att.bar.eastern <- mean(att.eastern)
    
    print(paste0("Eastern t-stat:",  att.bar.eastern , ", , t0:  ", t0))
    
    # bootstrap variance estimation
    df_mc_placebo_eastern <- widetoLong(Y= outcomes.placebo$M[!rownames(outcomes.placebo$M)%in%outcome.placebos$swiss,], 
                                        mask = (1-outcomes.placebo$mask[!rownames(outcomes.placebo$M)%in%outcome.placebos$swiss,]), 
                                        X = outcomes.placebo$X[!rownames(outcomes.placebo$M)%in%outcome.placebos$swiss,])
    df_mc_placebo_eastern$person_id <- as.numeric(factor(df_mc_placebo_eastern$person_id))
    att.eastern.boots <- clustered_bootstrap(current_data_realized_long=df_mc_placebo_eastern, estimator="mc_weights_covars", 
                                             N=nrow(outcomes.placebo$mask[!rownames(outcomes.placebo$M)%in%outcome.placebos$swiss,]), 
                                             T=ncol(outcomes.placebo$mask[!rownames(outcomes.placebo$M)%in%outcome.placebos$swiss,]), 
                                             B = 999, est_weights = TRUE, return_per_period=FALSE, return_replicates=TRUE)
    print(paste("Eastern variance:", round(var(att.eastern.boots),5)))
    
    boot.t.eastern.null <- att.eastern.boots - mean(att.eastern.boots,na.rm = TRUE) # center around zero
    
    boot.trajectory.eastern.pval <- (1+sum( abs(boot.t.eastern.null) > abs(att.bar.eastern)))/(999+1)
    
    print(paste0("Eastern p-val:", boot.trajectory.eastern.pval, ", t0: ", t0))
    
    # swiss
    
    att.swiss <- apply(mc.estimates.swiss.placebo$impact[which(rownames(mc.estimates.swiss.placebo$impact) %in% outcomes.placebo$swiss),]*outcomes.placebo$mask[which(rownames(outcomes.placebo$mask) %in% outcomes.placebo$swiss),],1,nzmean)[outcomes.placebo$swiss] 
    att.bar.swiss <- mean(att.swiss)
    
    print(paste0("swiss t-stat:",  att.bar.swiss , ", , t0:  ", t0))
    
    # bootstrap variance estimation
    df_mc_placebo_swiss <- widetoLong(Y= outcomes.placebo$M[which(rownames(outcomes.placebo$M) %in% outcomes.placebo$eastern),], 
                                        mask = (1-outcomes.placebo$mask[which(rownames(outcomes.placebo$M) %in% outcomes.placebo$eastern),]), 
                                        X = outcomes.placebo$X[which(rownames(outcomes.placebo$X) %in% outcomes.placebo$swiss),])
    df_mc_placebo_swiss$person_id <- as.numeric(factor(df_mc_placebo_swiss$person_id))
    att.swiss.boots <- clustered_bootstrap(current_data_realized_long=df_mc_placebo_swiss, estimator="mc_weights_covars", 
                                             N=nrow(outcomes.placebo$mask[which(rownames(outcomes.placebo$M) %in% outcomes.placebo$eastern),]), 
                                             T=ncol(outcomes.placebo$mask[which(rownames(outcomes.placebo$M) %in% outcomes.placebo$eastern),]), 
                                             B = 999, est_weights = TRUE, return_per_period=FALSE, return_replicates=TRUE)
    print(paste("swiss variance:", round(var(att.swiss.boots),5)))
    
    boot.t.swiss.null <- att.swiss.boots - mean(att.swiss.boots,na.rm = TRUE) # center around zero
    
    boot.trajectory.swiss.pval <- (1+sum( abs(boot.t.swiss.null) > abs(att.bar.swiss)))/(999+1)
    
    print(paste0("swiss p-val:", boot.trajectory.swiss.pval, ", t0: ", t0))
    
    print(paste0("placebo T0 ratio:," T0/ncol(outcomes.placebo$mask)))
    
  }
}