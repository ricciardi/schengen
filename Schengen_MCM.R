###################################
#  MC estimates with covariates   #
###################################

## Libraries
library(MCPanel)
library(Matrix)
library(data.table)
library(reshape)
library(reshape2)
library(matrixStats)
library(plyr)

source('utils.R')
source('MCEst.R')

## Set random seed
set.seed(10, "L'Ecuyer-CMRG")

outcome.vars <- c("CBWbord","CBWbordEMPL")

for(o in outcome.vars){
  print(paste0("Outcome: ", o))
  
  print(paste0("Covariates + FEs, outcome:",o))
  
  outcomes <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  # Get treatment effect estimates
    
  mc.estimates.eastern <- MCEst(outcomes, cluster='eastern', rev=TRUE, covars=TRUE)
  saveRDS(mc.estimates.eastern, paste0("results/mc-estimates-eastern-",o,".rds"))
  
  mc.estimates.swiss <- MCEst(outcomes, cluster='swiss', rev=TRUE, covars=TRUE)
  saveRDS(mc.estimates.swiss, paste0("results/mc-estimates-swiss-",o,".rds"))
  
  print(paste0("Rank of L (eastern): ", mc.estimates.eastern$rankL))
  print(paste0("Rank of L (swiss): ", mc.estimates.swiss$rankL))
  
  att.eastern <- apply(mc.estimates.eastern$impact[which(rownames(mc.estimates.eastern$impact) %in% outcomes$eastern),]*outcomes$mask[which(rownames(outcomes$mask) %in% outcomes$eastern),],1,nzmean)[outcomes$eastern] 
  att.bar.eastern <- mean(att.eastern)
  
  print(paste0("Eastern t-stat: Combined treatment effect (20051-20104): ", att.bar.eastern))
  
  att.swiss <- apply(mc.estimates.swiss$impact[which(rownames(mc.estimates.swiss$impact) %in% outcomes$swiss),]*outcomes$mask[which(rownames(outcomes$mask) %in% outcomes$swiss),],1,nzmean)[outcomes$swiss] 
  att.bar.swiss <- mean(att.swiss)
  
  print(paste0("swiss t-stat: Combined treatment effect (20051-20104): ", att.bar.swiss))
  
  # Eastern bootstrap variance estimation
  
  # bootstrap variance estimation
  df_mc_eastern <- widetoLong(Y= outcomes$M[!rownames(outcomes$M)%in%outcomes$swiss,], 
                      mask = (1-outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]), 
                      X = outcomes$X[!rownames(outcomes$M)%in%outcomes$swiss,])
  df_mc_eastern$person_id <- as.numeric(factor(df_mc_eastern$person_id))
  att.eastern.boot.per.period <- clustered_bootstrap(current_data_realized_long=df_mc_eastern, 
                                                     estimator="mc_weights_covars", 
                                                     N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]), 
                                                     T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]),
                                                     B = 999, 
                                                     est_weights = TRUE, return_per_period=TRUE, return_replicates=FALSE)
  saveRDS(att.eastern.boot.per.period, paste0("results/boot-eastern-",o,".rds")) 
  
  print("Eastern CI: Combined treatment effect (20051-20104):")
  print(unlist(boot_CI(att.bar.eastern, att.eastern.boot.per.period$att.bar.boot.var)))
  
  att.eastern.boots <- clustered_bootstrap(current_data_realized_long=df_mc_eastern, estimator="mc_weights_covars", 
                                           N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]), 
                                           T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]), B = 999, est_weights = TRUE, return_per_period=FALSE, return_replicates=TRUE)
  
  boot.t.eastern.null <- att.eastern.boots - mean(att.eastern.boots,na.rm = TRUE) # center around zero
  
  boot.trajectory.eastern.pval <- (1+sum( abs(boot.t.eastern.null) > abs(att.bar.eastern)))/(999+1)
  
  print(paste0("Eastern p-val: Combined treatment effect (20051-20104): ", boot.trajectory.eastern.pval))
  
  print(paste0("Eastern effect share: Combined treatment effect (20051-20104): ", (att.bar.eastern)/mean(outcomes$M[rownames(outcomes$M)%in%outcomes$eastern,][,1:(which(colnames(outcomes$mask)==20111)-1)])))
  
  saveRDS(att.eastern.boots, paste0("results/boot-trajectory-eastern-",o,".rds")) 
  
  # Swiss variance estimation
  
  # bootstrap variance estimation
  df_mc_swiss <- widetoLong(Y= outcomes$M[!rownames(outcomes$M)%in%outcomes$eastern,], 
                              mask = (1-outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]), 
                              X = outcomes$X[!rownames(outcomes$M)%in%outcomes$eastern,])
  df_mc_swiss$person_id <- as.numeric(factor(df_mc_swiss$person_id))
  att.swiss.boot.per.period <- clustered_bootstrap(current_data_realized_long=df_mc_swiss, estimator="mc_weights_covars", 
                                                     N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]), 
                                                     T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]),
                                                     B = 999, est_weights = TRUE, return_per_period=TRUE, return_replicates=FALSE)
  saveRDS(att.swiss.boot.per.period, paste0("results/boot-swiss-",o,".rds")) 
  
  print("Swiss CI: Combined treatment effect (20051-20104):")
  print(unlist(boot_CI(att.bar.swiss, att.swiss.boot.per.period$att.bar.boot.var)))
  
  att.swiss.boots <- clustered_bootstrap(current_data_realized_long=df_mc_swiss, estimator="mc_weights_covars", 
                                           N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]), 
                                           T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]), B = 999, est_weights = TRUE, return_per_period=FALSE, return_replicates=TRUE)
  
  boot.t.swiss.null <- att.swiss.boots - mean(att.swiss.boots,na.rm = TRUE) # center around zero
  
  boot.trajectory.swiss.pval <- (1+sum( abs(boot.t.swiss.null) > abs(att.bar.swiss)))/(999+1)
  
  print(paste0("Swiss p-val: Combined treatment effect (20051-20104): ", boot.trajectory.swiss.pval))
  
  print(paste0("Swiss effect share: Combined treatment effect (20051-20104): ", (att.bar.swiss)/mean(outcomes$M[rownames(outcomes$M)%in%outcomes$swiss,][,1:(which(colnames(outcomes$mask)==20111)-1)])))
  
  saveRDS(att.eastern.boots, paste0("results/boot-trajectory-swiss-",o,".rds")) 
  
  ## Eastern: T0' = 20081
  
  outcomes.20081 <- outcomes
  outcomes.20081$mask[rownames(outcomes.20081$mask) %in% outcomes.20081$eastern,][,which(colnames(outcomes.20081$mask)=="20081"):ncol(outcomes.20081$mask)] <- 0
  
  mc.estimates.1 <- MCEst(outcomes.20081, cluster='eastern', rev=TRUE, covars=TRUE) 

  att.eastern.1 <- apply(mc.estimates.1$impact[which(rownames(mc.estimates.1$impact) %in% outcomes.20081$eastern),]*outcomes.20081$mask[which(rownames(outcomes.20081$mask) %in% outcomes.20081$eastern),],1,nzmean)[outcomes.20081$eastern] 
  att.bar.eastern.1 <- mean(att.eastern.1)
  
  print(paste0("Eastern t-stat: partial treatment effect (20051-20074): ", att.bar.eastern.1))
  
  # Eastern (1) bootstrap variance estimation
  
  # bootstrap variance estimation
  df_mc_eastern_1 <- widetoLong(Y= outcomes.20081$M[!rownames(outcomes.20081$M)%in%outcomes.20081$swiss,], 
                              mask = (1-outcomes.20081$mask[!rownames(outcomes.20081$M)%in%outcomes.20081$swiss,]), 
                              X = outcomes.20081$X[!rownames(outcomes.20081$M)%in%outcomes.20081$swiss,])
  df_mc_eastern_1$person_id <- as.numeric(factor(df_mc_eastern_1$person_id))
  att.eastern.1.boot.per.period <- clustered_bootstrap(current_data_realized_long=df_mc_eastern_1, estimator="mc_weights_covars", 
                                                     N=nrow(outcomes.20081$mask[!rownames(outcomes.20081$M)%in%outcomes.20081$swiss,]), 
                                                     T=ncol(outcomes.20081$mask[!rownames(outcomes.20081$M)%in%outcomes.20081$swiss,]),
                                                     B = 999, est_weights = TRUE, return_per_period=TRUE, return_replicates=FALSE)

  print("Eastern CI: Combined treatment effect (20051-20074):")
  print(unlist(boot_CI(att.bar.eastern.1, att.eastern.1.boot.per.period$att.bar.boot.var)))
  
  att.eastern.1.boots <- clustered_bootstrap(current_data_realized_long=df_mc_eastern_1, estimator="mc_weights_covars", 
                                           N=nrow(outcomes.20081$mask[!rownames(outcomes.20081$M)%in%outcomes.20081$swiss,]), 
                                           T=ncol(outcomes.20081$mask[!rownames(outcomes.20081$M)%in%outcomes.20081$swiss,]), B = 999, est_weights = TRUE, return_per_period=FALSE, return_replicates=TRUE)
  
  boot.t.eastern.1.null <- att.eastern.1.boots - mean(att.eastern.1.boots,na.rm = TRUE) # center around zero
  
  boot.trajectory.eastern.1.pval <- (1+sum( abs(boot.t.eastern.1.null) > abs(att.bar.eastern.1)))/(999+1)
  
  print(paste0("Eastern p-val: Combined treatment effect (20051-20074): ", boot.trajectory.eastern.1.pval))
  
  ## Eastern: T0' = 20111
  
  outcomes.20111 <- outcomes
  outcomes.20111$mask[rownames(outcomes.20111$mask) %in% outcomes.20111$eastern,][,which(colnames(outcomes.20111$mask)=="20111"):ncol(outcomes.20111$mask)] <- 0
  
  mc.estimates.2 <- MCEst(outcomes.20111, cluster='eastern', rev=TRUE, covars=TRUE) 
  
  att.eastern.2 <- apply(mc.estimates.2$impact[which(rownames(mc.estimates.2$impact) %in% outcomes.20111$eastern),]*outcomes.20111$mask[which(rownames(outcomes.20111$mask) %in% outcomes.20111$eastern),],1,nzmean)[outcomes.20111$eastern] 
  att.bar.eastern.2 <- mean(att.eastern.2)
  
  print(paste0("Eastern t-stat: partial treatment effect (20081-20104): ", att.bar.eastern.2))
  
  # Eastern (2) bootstrap variance estimation
  
  # bootstrap variance estimation
  df_mc_eastern_2 <- widetoLong(Y= outcomes.20111$M[!rownames(outcomes.20111$M)%in%outcomes.20111$swiss,], 
                                mask = (1-outcomes.20111$mask[!rownames(outcomes.20111$M)%in%outcomes.20111$swiss,]), 
                                X = outcomes.20111$X[!rownames(outcomes.20111$X)%in%outcomes.20111$swiss,])
  df_mc_eastern_2$person_id <- as.numeric(factor(df_mc_eastern_2$person_id))
  att.eastern.2.boot.per.period <- clustered_bootstrap(current_data_realized_long=df_mc_eastern_2, estimator="mc_weights_covars", 
                                                       N=nrow(outcomes.20111$mask[!rownames(outcomes.20111$M)%in%outcomes.20111$swiss,]), 
                                                       T=ncol(outcomes.20111$mask[!rownames(outcomes.20111$M)%in%outcomes.20111$swiss,]),
                                                       B = 999, est_weights = TRUE, return_per_period=TRUE, return_replicates=FALSE)
  
  print("Eastern CI: Combined treatment effect (20081-20104):")
  print(unlist(boot_CI(att.bar.eastern.2, att.eastern.2.boot.per.period$att.bar.boot.var)))
  
  att.eastern.2.boots <- clustered_bootstrap(current_data_realized_long=df_mc_eastern_2, estimator="mc_weights_covars", 
                                             N=nrow(outcomes.20111$mask[!rownames(outcomes.20111$M)%in%outcomes.20111$swiss,]), 
                                             T=ncol(outcomes.20111$mask[!rownames(outcomes.20111$M)%in%outcomes.20111$swiss,]), B = 999, est_weights = TRUE, return_per_period=FALSE, return_replicates=TRUE)
  
  boot.t.eastern.2.null <- att.eastern.2.boots - mean(att.eastern.2.boots,na.rm = TRUE) # center around zero
  
  boot.trajectory.eastern.2.pval <- (1+sum( abs(boot.t.eastern.2.null) > abs(att.bar.eastern.2)))/(999+1)
  
  print(paste0("Eastern p-val: Combined treatment effect (20081-20104): ", boot.trajectory.eastern.2.pval))
  
  ###
  ## swiss: T0' = 20071
  
  outcomes.20071 <- outcomes
  outcomes.20071$mask[rownames(outcomes.20071$mask) %in% outcomes.20071$swiss,][,which(colnames(outcomes.20071$mask)=="20071"):ncol(outcomes.20071$mask)] <- 0
  
  mc.estimates.swiss.1 <- MCEst(outcomes.20071, cluster='swiss', rev=TRUE, covars=TRUE) 
  
  att.swiss.1 <- apply(mc.estimates.swiss.1$impact[which(rownames(mc.estimates.swiss.1$impact) %in% outcomes.20071$swiss),]*outcomes.20071$mask[which(rownames(outcomes.20071$mask) %in% outcomes.20071$swiss),],1,nzmean)[outcomes.20071$swiss] 
  att.bar.swiss.1 <- mean(att.swiss.1)
  
  print(paste0("swiss t-stat: partial treatment effect (20051-20072): ", att.bar.swiss.1))
  
  # swiss (1) bootstrap variance estimation
  
  # bootstrap variance estimation
  df_mc_swiss_1 <- widetoLong(Y= outcomes.20071$M[!rownames(outcomes.20071$M)%in%outcomes.20071$eastern,], 
                                mask = (1-outcomes.20071$mask[!rownames(outcomes.20071$M)%in%outcomes.20071$eastern,]), 
                                X = outcomes.20071$X[!rownames(outcomes.20071$M)%in%outcomes.20071$eastern,])
  df_mc_swiss_1$person_id <- as.numeric(factor(df_mc_swiss_1$person_id))
  att.swiss.1.boot.per.period <- clustered_bootstrap(current_data_realized_long=df_mc_swiss_1, estimator="mc_weights_covars", 
                                                       N=nrow(outcomes.20071$mask[!rownames(outcomes.20071$M)%in%outcomes.20071$eastern,]), 
                                                       T=ncol(outcomes.20071$mask[!rownames(outcomes.20071$M)%in%outcomes.20071$eastern,]),
                                                       B = 999, est_weights = TRUE, return_per_period=TRUE, return_replicates=FALSE)
  
  print("swiss CI: Combined treatment effect (20051-20072):")
  print(unlist(boot_CI(att.bar.swiss.1, att.swiss.1.boot.per.period$att.bar.boot.var)))
  
  att.swiss.1.boots <- clustered_bootstrap(current_data_realized_long=df_mc_swiss_1, estimator="mc_weights_covars", 
                                             N=nrow(outcomes.20071$mask[!rownames(outcomes.20071$M)%in%outcomes.20071$eastern,]), 
                                             T=ncol(outcomes.20071$mask[!rownames(outcomes.20071$M)%in%outcomes.20071$eastern,]), B = 999, est_weights = TRUE, return_per_period=FALSE, return_replicates=TRUE)
  
  boot.t.swiss.1.null <- att.swiss.1.boots - mean(att.swiss.1.boots,na.rm = TRUE) # center around zero
  
  boot.trajectory.swiss.1.pval <- (1+sum( abs(boot.t.swiss.1.null) > abs(att.bar.swiss.1)))/(999+1)
  
  print(paste0("swiss p-val: Combined treatment effect (20051-20072): ", boot.trajectory.swiss.1.pval))
  
  ## swiss: T0' = 20091
  
  outcomes.20091 <- outcomes
  outcomes.20091$mask[rownames(outcomes.20091$mask) %in% outcomes.20091$swiss,][,which(colnames(outcomes.20091$mask)=="20091"):ncol(outcomes.20091$mask)] <- 0
  
  mc.estimates.swiss.2 <- MCEst(outcomes.20091, cluster='swiss', rev=TRUE, covars=TRUE) 
  
  att.swiss.2 <- apply(mc.estimates.swiss.2$impact[which(rownames(mc.estimates.swiss.2$impact) %in% outcomes.20091$swiss),]*outcomes.20091$mask[which(rownames(outcomes.20091$mask) %in% outcomes.20091$swiss),],1,nzmean)[outcomes.20091$swiss] 
  att.bar.swiss.2 <- mean(att.swiss.2)
  
  print(paste0("swiss t-stat: partial treatment effect (20073-20084): ", att.bar.swiss.2))
  
  # swiss (2) bootstrap variance estimation
  
  # bootstrap variance estimation
  df_mc_swiss_2 <- widetoLong(Y= outcomes.20091$M[!rownames(outcomes.20091$M)%in%outcomes.20091$eastern,], 
                                mask = (1-outcomes.20091$mask[!rownames(outcomes.20091$M)%in%outcomes.20091$eastern,]), 
                                X = outcomes.20091$X[!rownames(outcomes.20091$X)%in%outcomes.20091$eastern,])
  df_mc_swiss_2$person_id <- as.numeric(factor(df_mc_swiss_2$person_id))
  att.swiss.2.boot.per.period <- clustered_bootstrap(current_data_realized_long=df_mc_swiss_2, estimator="mc_weights_covars", 
                                                       N=nrow(outcomes.20091$mask[!rownames(outcomes.20091$M)%in%outcomes.20091$eastern,]), 
                                                       T=ncol(outcomes.20091$mask[!rownames(outcomes.20091$M)%in%outcomes.20091$eastern,]),
                                                       B = 999, est_weights = TRUE, return_per_period=TRUE, return_replicates=FALSE)
  
  print("swiss CI: Combined treatment effect (20073-20084):")
  print(unlist(boot_CI(att.bar.swiss.2, att.swiss.2.boot.per.period$att.bar.boot.var)))
  
  att.swiss.2.boots <- clustered_bootstrap(current_data_realized_long=df_mc_swiss_2, estimator="mc_weights_covars", 
                                             N=nrow(outcomes.20091$mask[!rownames(outcomes.20091$M)%in%outcomes.20091$eastern,]), 
                                             T=ncol(outcomes.20091$mask[!rownames(outcomes.20091$M)%in%outcomes.20091$eastern,]), B = 999, est_weights = TRUE, return_per_period=FALSE, return_replicates=TRUE)
  
  boot.t.swiss.2.null <- att.swiss.2.boots - mean(att.swiss.2.boots,na.rm = TRUE) # center around zero
  
  boot.trajectory.swiss.2.pval <- (1+sum( abs(boot.t.swiss.2.null) > abs(att.bar.swiss.2)))/(999+1)
  
  print(paste0("swiss p-val: Combined treatment effect (20073-20084): ", boot.trajectory.swiss.2.pval))
  }