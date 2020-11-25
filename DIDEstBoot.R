DIDEstBoot <- function(tseries,mask,W,X=NULL,X.hat=NULL, eastern=NULL, swiss=NULL, est_eastern=FALSE, est_swiss=FALSE, rev=TRUE, estimator="DID") {
  
  if(est_eastern){
    Y <- t(tseries)[!rownames(t(tseries))%in%swiss,] # NxT
  }else{
    Y <- t(tseries)[!rownames(t(tseries))%in%eastern,] # NxT
  }
  
  treat <- mask # NxT masked matrix 
  treat <- treat[rownames(treat) %in% row.names(Y),]
  
  N <- nrow(treat)
  T <- ncol(treat)
  
  treat_mat <- 1-treat
  treat_mat_NA <- treat_mat
  treat_mat_NA[treat_mat_NA==0] <- NA # zeros are NA
  
  Y_obs <- Y * treat_mat
  Y_obs_NA <- Y * treat_mat_NA
  
  if(estimator=="DID"){
    
    est_model <- t(DID(t(Y_obs), t(treat_mat)))
    
  }
  
  if(estimator=="ADH"){
    
    est_model <- adh_mp_rows(Y_obs, treat_mat, rel_tol = 0.001)
    
  }
  
  if(estimator=="ENT"){
    
    est_model <- t(en_mp_rows(t(Y_obs), t(treat_mat), num_folds = 3))
  }
  
  if(estimator=="NNMF"){
    nnmf_fit <- nnmf(Y_obs_NA,1,verbose=0)
    est_model <- Y_obs + with(nnmf_fit, W %*% H)*treat # impute only missing
  }
  
  if(rev){
    impact <- (est_model-Y)
    return(impact)
  }else{
    impact <- (Y-est_model)
    return(impact)
  }
}