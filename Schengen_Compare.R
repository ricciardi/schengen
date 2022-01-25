###################################
#  DID, IFE, SCM estimates #
###################################

library(MCPanel)
library(Matrix)
library(data.table)
library(reshape)
library(reshape2)
library(matrixStats)
library(plyr)
library(emfactor)

source('utils.R')
source('MCEst.R')
source('IFE.R')

## Set random seed
set.seed(10, "L'Ecuyer-CMRG")

outcome.vars <- c("CBWbord","CBWbordEMPL")

for(o in outcome.vars){
  print(paste0("Outcome: ", o))
  
  outcomes <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  # ADH Estimation 
  
  # Get treatment effect estimates
  
  adh.estimates.eastern <- MCEst(outcomes, cluster='eastern', rev=TRUE, covars=FALSE, ADH=TRUE)
  
  adh.estimates.swiss <- MCEst(outcomes, cluster='swiss', rev=TRUE, covars=FALSE, ADH=TRUE)
  
  adh.att.eastern <- apply(adh.estimates.eastern$impact[which(rownames(adh.estimates.eastern$impact) %in% outcomes$eastern),]*outcomes$mask[which(rownames(outcomes$mask) %in% outcomes$eastern),],1,nzmean)[outcomes$eastern] 
  adh.adh.att.bar.eastern <- mean(adh.att.eastern)
  
  print(paste0("ADH Eastern t-stat: Combined treatment effect (20051-20104): ", adh.adh.att.bar.eastern))
  
  adh.att.swiss <- apply(adh.estimates.swiss$impact[which(rownames(adh.estimates.swiss$impact) %in% outcomes$swiss),]*outcomes$mask[which(rownames(outcomes$mask) %in% outcomes$swiss),],1,nzmean)[outcomes$swiss] 
  adh.adh.att.bar.swiss <- mean(adh.att.swiss)
  
  print(paste0("ADH swiss t-stat: Combined treatment effect (20051-20104): ", adh.adh.att.bar.swiss))
  
  # Eastern bootstrap variance estimation
  
  # bootstrap variance estimation
  df_eastern <- widetoLong(Y= outcomes$M[!rownames(outcomes$M)%in%outcomes$swiss,], 
                              mask = (1-outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]), 
                              X = outcomes$X[!rownames(outcomes$M)%in%outcomes$swiss,])
  df_eastern$person_id <- as.numeric(factor(df_eastern$person_id))
  adh.att.eastern.boot.per.period <- clustered_bootstrap(current_data_realized_long=df_eastern, 
                                                     estimator="ADH", 
                                                     N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]), 
                                                     T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]),
                                                     B = 999, 
                                                     est_weights = TRUE, return_per_period=TRUE, return_replicates=FALSE)
  
  print("ADH Eastern CI: Combined treatment effect (20051-20104):")
  print(unlist(boot_CI(adh.adh.att.bar.eastern, adh.att.eastern.boot.per.period$att.bar.boot.var)))
  
  adh.att.eastern.boots <- clustered_bootstrap(current_data_realized_long=df_eastern, estimator="ADH", 
                                           N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]), 
                                           T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]), B = 999, est_weights = TRUE, return_per_period=FALSE, return_replicates=TRUE)
  
  adh.boot.t.eastern.null <- adh.att.eastern.boots - mean(adh.att.eastern.boots,na.rm = TRUE) # center around zero
  
  boot.trajectory.eastern.pval <- (1+sum( abs(adh.boot.t.eastern.null) > abs(adh.att.bar.eastern)))/(999+1)
  
  print(paste0("ADH Eastern p-val: Combined treatment effect (20051-20104): ", boot.trajectory.eastern.pval))
  
  print(paste0("ADH Eastern effect share: Combined treatment effect (20051-20104): ", (adh.att.bar.eastern)/mean(outcomes$M[rownames(outcomes$M)%in%outcomes$eastern,][,1:(which(colnames(outcomes$mask)==20111)-1)])))
  
  # Swiss variance estimation
  
  # bootstrap variance estimation
  df_swiss <- widetoLong(Y= outcomes$M[!rownames(outcomes$M)%in%outcomes$eastern,], 
                            mask = (1-outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]), 
                            X = outcomes$X[!rownames(outcomes$M)%in%outcomes$eastern,])
  df_swiss$person_id <- as.numeric(factor(df_swiss$person_id))
  adh.att.swiss.boot.per.period <- clustered_bootstrap(current_data_realized_long=df_swiss, estimator="ADH", 
                                                   N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]), 
                                                   T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]),
                                                   B = 999, est_weights = TRUE, return_per_period=TRUE, return_replicates=FALSE)
  
  print("ADH Swiss CI: Combined treatment effect (20051-20104):")
  print(unlist(boot_CI(adh.att.bar.swiss, adh.att.swiss.boot.per.period$att.bar.boot.var)))
  
  adh.att.swiss.boots <- clustered_bootstrap(current_data_realized_long=df_swiss, estimator="ADH", 
                                         N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]), 
                                         T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]), B = 999, est_weights = TRUE, return_per_period=FALSE, return_replicates=TRUE)
  
  adh.boot.t.swiss.null <- adh.att.swiss.boots - mean(adh.att.swiss.boots,na.rm = TRUE) # center around zero
  
  boot.trajectory.swiss.pval <- (1+sum( abs(adh.boot.t.swiss.null) > abs(adh.att.bar.swiss)))/(999+1)
  
  print(paste0("ADH Swiss p-val: Combined treatment effect (20051-20104): ", boot.trajectory.swiss.pval))
  
  print(paste0("ADH Swiss effect share: Combined treatment effect (20051-20104): ", (adh.att.bar.swiss)/mean(outcomes$M[rownames(outcomes$M)%in%outcomes$swiss,][,1:(which(colnames(outcomes$mask)==20111)-1)])))
  
  # DID Estimation 
  
  # Get treatment effect estimates
  
  did.estimates.eastern <- MCEst(outcomes, cluster='eastern', rev=TRUE, covars=FALSE, DID=TRUE)
  
  did.estimates.swiss <- MCEst(outcomes, cluster='swiss', rev=TRUE, covars=FALSE, DID=TRUE)
  
  did.att.eastern <- apply(did.estimates.eastern$impact[which(rownames(did.estimates.eastern$impact) %in% outcomes$eastern),]*outcomes$mask[which(rownames(outcomes$mask) %in% outcomes$eastern),],1,nzmean)[outcomes$eastern] 
  did.att.bar.eastern <- mean(did.att.eastern)
  
  print(paste0("DID Eastern t-stat: Combined treatment effect (20051-20104): ", did.att.bar.eastern))
  
  did.att.swiss <- apply(did.estimates.swiss$impact[which(rownames(did.estimates.swiss$impact) %in% outcomes$swiss),]*outcomes$mask[which(rownames(outcomes$mask) %in% outcomes$swiss),],1,nzmean)[outcomes$swiss] 
  did.att.bar.swiss <- mean(did.att.swiss)
  
  print(paste0("DID swiss t-stat: Combined treatment effect (20051-20104): ", did.att.bar.swiss))
  
  # Eastern bootstrap variance estimation
  
  # bootstrap variance estimation
  did.att.eastern.boot.per.period <- clustered_bootstrap(current_data_realized_long=df_eastern, 
                                                         estimator="DID", 
                                                         N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]), 
                                                         T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]),
                                                         B = 999, 
                                                         est_weights = TRUE, return_per_period=TRUE, return_replicates=FALSE)
  
  print("DID Eastern CI: Combined treatment effect (20051-20104):")
  print(unlist(boot_CI(did.att.bar.eastern, did.att.eastern.boot.per.period$att.bar.boot.var)))
  
  did.att.eastern.boots <- clustered_bootstrap(current_data_realized_long=df_eastern, estimator="DID", 
                                               N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]), 
                                               T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]), B = 999, est_weights = TRUE, return_per_period=FALSE, return_replicates=TRUE)
  
  did.boot.t.eastern.null <- did.att.eastern.boots - mean(did.att.eastern.boots,na.rm = TRUE) # center around zero
  
  boot.trajectory.eastern.pval <- (1+sum( abs(did.boot.t.eastern.null) > abs(did.att.bar.eastern)))/(999+1)
  
  print(paste0("DID Eastern p-val: Combined treatment effect (20051-20104): ", boot.trajectory.eastern.pval))
  
  print(paste0("DID Eastern effect share: Combined treatment effect (20051-20104): ", (did.att.bar.eastern)/mean(outcomes$M[rownames(outcomes$M)%in%outcomes$eastern,][,1:(which(colnames(outcomes$mask)==20111)-1)])))
  
  # Swiss variance estimation
  
  # bootstrap variance estimation
  did.att.swiss.boot.per.period <- clustered_bootstrap(current_data_realized_long=df_swiss, estimator="DID", 
                                                       N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]), 
                                                       T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]),
                                                       B = 999, est_weights = TRUE, return_per_period=TRUE, return_replicates=FALSE)
  
  print("DID Swiss CI: Combined treatment effect (20051-20104):")
  print(unlist(boot_CI(did.att.bar.swiss, did.att.swiss.boot.per.period$att.bar.boot.var)))
  
  did.att.swiss.boots <- clustered_bootstrap(current_data_realized_long=df_swiss, estimator="DID", 
                                             N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]), 
                                             T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]), B = 999, est_weights = TRUE, return_per_period=FALSE, return_replicates=TRUE)
  
  did.boot.t.swiss.null <- did.att.swiss.boots - mean(did.att.swiss.boots,na.rm = TRUE) # center around zero
  
  boot.trajectory.swiss.pval <- (1+sum( abs(did.boot.t.swiss.null) > abs(did.att.bar.swiss)))/(999+1)
  
  print(paste0("DID Swiss p-val: Combined treatment effect (20051-20104): ", boot.trajectory.swiss.pval))
  
  print(paste0("DID Swiss effect share: Combined treatment effect (20051-20104): ", (did.att.bar.swiss)/mean(outcomes$M[rownames(outcomes$M)%in%outcomes$swiss,][,1:(which(colnames(outcomes$mask)==20111)-1)])))
  
  # IFE Estimation 
  
  # Get treatment effect estimates
  
  ife.estimates.eastern <- MCEst(outcomes, cluster='eastern', rev=TRUE, covars=FALSE, IFE=TRUE)
  
  ife.estimates.swiss <- MCEst(outcomes, cluster='swiss', rev=TRUE, covars=FALSE, IFE=TRUE)
  
  ife.att.eastern <- apply(ife.estimates.eastern$impact[which(rownames(ife.estimates.eastern$impact) %in% outcomes$eastern),]*outcomes$mask[which(rownames(outcomes$mask) %in% outcomes$eastern),],1,nzmean)[outcomes$eastern] 
  ife.att.bar.eastern <- mean(ife.att.eastern)
  
  print(paste0("IFE Eastern t-stat: Combined treatment effect (20051-20104): ", ife.att.bar.eastern))
  
  ife.att.swiss <- apply(ife.estimates.swiss$impact[which(rownames(ife.estimates.swiss$impact) %in% outcomes$swiss),]*outcomes$mask[which(rownames(outcomes$mask) %in% outcomes$swiss),],1,nzmean)[outcomes$swiss] 
  ife.att.bar.swiss <- mean(ife.att.swiss)
  
  print(paste0("IFE swiss t-stat: Combined treatment effect (20051-20104): ", ife.att.bar.swiss))
  
  # Eastern bootstrap variance estimation
  
  # bootstrap variance estimation
  ife.att.eastern.boot.per.period <- clustered_bootstrap(current_data_realized_long=df_eastern, 
                                                         estimator="IFE", 
                                                         N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]), 
                                                         T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]),
                                                         B = 999, 
                                                         est_weights = TRUE, return_per_period=TRUE, return_replicates=FALSE)
  
  print("IFE Eastern CI: Combined treatment effect (20051-20104):")
  print(unlist(boot_CI(ife.att.bar.eastern, ife.att.eastern.boot.per.period$att.bar.boot.var)))
  
  ife.att.eastern.boots <- clustered_bootstrap(current_data_realized_long=df_eastern, estimator="IFE", 
                                               N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]), 
                                               T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$swiss,]), B = 999, est_weights = TRUE, return_per_period=FALSE, return_replicates=TRUE)
  
  ife.boot.t.eastern.null <- ife.att.eastern.boots - mean(ife.att.eastern.boots,na.rm = TRUE) # center around zero
  
  boot.trajectory.eastern.pval <- (1+sum( abs(ife.boot.t.eastern.null) > abs(ife.att.bar.eastern)))/(999+1)
  
  print(paste0("IFE Eastern p-val: Combined treatment effect (20051-20104): ", boot.trajectory.eastern.pval))
  
  print(paste0("IFE Eastern effect share: Combined treatment effect (20051-20104): ", (ife.att.bar.eastern)/mean(outcomes$M[rownames(outcomes$M)%in%outcomes$eastern,][,1:(which(colnames(outcomes$mask)==20111)-1)])))
  
  # Swiss variance estimation
  
  # bootstrap variance estimation
  ife.att.swiss.boot.per.period <- clustered_bootstrap(current_data_realized_long=df_swiss, estimator="IFE", 
                                                       N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]), 
                                                       T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]),
                                                       B = 999, est_weights = TRUE, return_per_period=TRUE, return_replicates=FALSE)
  
  print("IFE Swiss CI: Combined treatment effect (20051-20104):")
  print(unlist(boot_CI(ife.att.bar.swiss, ife.att.swiss.boot.per.period$att.bar.boot.var)))
  
  ife.att.swiss.boots <- clustered_bootstrap(current_data_realized_long=df_swiss, estimator="IFE", 
                                             N=nrow(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]), 
                                             T=ncol(outcomes$mask[!rownames(outcomes$M)%in%outcomes$eastern,]), B = 999, est_weights = TRUE, return_per_period=FALSE, return_replicates=TRUE)
  
  ife.boot.t.swiss.null <- ife.att.swiss.boots - mean(ife.att.swiss.boots,na.rm = TRUE) # center around zero
  
  boot.trajectory.swiss.pval <- (1+sum( abs(ife.boot.t.swiss.null) > abs(ife.att.bar.swiss)))/(999+1)
  
  print(paste0("IFE Swiss p-val: Combined treatment effect (20051-20104): ", boot.trajectory.swiss.pval))
  
  print(paste0("IFE Swiss effect share: Combined treatment effect (20051-20104): ", (ife.att.bar.swiss)/mean(outcomes$M[rownames(outcomes$M)%in%outcomes$swiss,][,1:(which(colnames(outcomes$mask)==20111)-1)])))
}