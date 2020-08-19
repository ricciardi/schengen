DIDEstBoot <- function(tseries,mask,W,X=NULL,X.hat=NULL, t0=NULL, eastern=NULL, swiss=NULL, control=NULL, covars=TRUE, rev=TRUE, fe=TRUE) {
  
  Y <- t(tseries) # NxT 
  
  treat <- mask # NxT masked matrix 
  
  N <- nrow(treat)
  T <- ncol(treat)
  
  treat_mat <- 1-treat
  
  Y_obs <- Y * treat_mat
  
  ## ------
  ##  DID
  ## ------
  
  est_model_DID <- t(DID(t(Y_obs), t(1-outcomes.cbw$mask)))
  
  if(rev){
    impact <- (est_model_DID-Y)
    
    if(!is.null(eastern)){
      att <- as.matrix(colMeans(impact[,1:(t0-1)][rownames(impact) %in% eastern,])) # get mean pre-period impact on treated
      return(att)
    }
    if(!is.null(swiss)){
      att <- as.matrix(colMeans(impact[,1:(t0-1)][rownames(impact) %in% swiss,]))
      return(att)
    } 
    if(!is.null(control)){
      att<- as.matrix(colMeans(impact[,1:(t0-1)][rownames(impact) %in% control,]))
      return(att)
    } else{
      return(impact)
    }
  }else{
    impact <- (Y-est_model_DID)
    if(!is.null(eastern)){
      att <- as.matrix(colMeans(impact[,t0:ncol(impact)][rownames(impact) %in% eastern,])) # get mean post-period impact on treated
      return(att)
    }
    if(!is.null(swiss)){
      as.matrix(colMeans(impact[,t0:ncol(impact)][rownames(impact) %in% swiss,]))
      return(att) 
    } 
    if(!is.null(control)){
      att <- as.matrix(colMeans(impact[,t0:ncol(impact)][rownames(impact) %in% control,]))
      return(att)
    }else{
      return(impact)
    }
  }
}