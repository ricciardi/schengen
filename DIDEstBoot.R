DIDEstBoot <- function(tseries,mask,W,X=NULL,X.hat=NULL, t0=NULL, eastern=NULL, swiss=NULL, rev=TRUE, fe=TRUE, estimator="DID") {
  
  Y <- t(tseries) # NxT 
  
  treat <- mask # NxT masked matrix 
  
  N <- nrow(treat)
  T <- ncol(treat)
  
  treat_mat <- 1-treat
  
  Y_obs <- Y * treat_mat
  
  if(estimator=="DID"){

    est_model <- t(DID(t(Y_obs), t(1-outcomes.cbw$mask)))
    
  }
  
  if(estimator=="ADH"){
    
    est_model <- adh_mp_rows(Y_obs, 1-outcomes.cbw$mask, rel_tol = 0.001)
    
  }
  
  if(estimator=="ENT"){
    
    est_model <- t(en_mp_rows(t(Y_obs), t(1-outcomes.cbw$mask), num_folds = 3))
  }
  
  if(rev){
    impact <- (est_model-Y)
    
    if(!is.null(eastern)){
      att <- as.matrix(colMeans(impact[,1:(t0-1)][rownames(impact) %in% eastern,])) # get mean pre-period impact on treated
      return(att)
    }
    if(!is.null(swiss)){
      att <- as.matrix(colMeans(impact[,1:(t0-1)][rownames(impact) %in% swiss,]))
      return(att)
    } else{
      return(impact)
    }
  }else{
    impact <- (Y-est_model)
    if(!is.null(eastern)){
      att <- as.matrix(colMeans(impact[,t0:ncol(impact)][rownames(impact) %in% eastern,])) # get mean post-period impact on treated
      return(att)
    }
    if(!is.null(swiss)){
      as.matrix(colMeans(impact[,t0:ncol(impact)][rownames(impact) %in% swiss,]))
      return(att) 
    } else{
      return(impact)
    }
  }
}