MCEstBoot <- function(tseries,mask,W,X=NULL,X.hat=NULL, t0=NULL, z.cbw.eastern, z.cbw.swiss, eastern=NULL, swiss=NULL, est_eastern=FALSE, est_swiss=FALSE, covars=TRUE, rev=TRUE, best_L, best_B) {
  
  Y <- t(tseries) # NxT 
  
  treat <- mask # NxT masked matrix 
  
  N <- nrow(treat)
  T <- ncol(treat)
  
  treat_mat <- 1-treat
  
  Y_obs <- Y * treat_mat
  
  W <- W
  W <- W[rownames(W) %in% row.names(Y),]
  W <- W[row.names(Y),]  # reorder
  W[W<=0] <- min(W[W>0]) # set floor
  W[W>=1] <- max(W[W<1]) # set ceiling
  
  weights <- matrix(NA, nrow=nrow(W), ncol=ncol(W), dimnames = list(rownames(W), colnames(W)))
  weights <- (1-treat) + (treat)*((1-W)/(W))  
  weights[rownames(weights) %in% eastern,] <- weights[rownames(weights) %in% eastern,] %*%diag(z.cbw.eastern)
  weights[rownames(weights) %in% swiss,] <- weights[rownames(weights) %in% swiss,] %*%diag(z.cbw.swiss)
  
  if(covars){
    
    ## ------
    ## MC-NNM-W
    ## ------
    
    est_model_MCPanel_w <- mcnnm_wc(M = Y_obs, C = X, mask = treat_mat, W = weights, to_normalize = 1, to_estimate_u = 1, to_estimate_v = 1, lambda_L=best_L, lambda_B = best_B, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]] 
      
    est_model_MCPanel_w$Mhat <- est_model_MCPanel_w$L + X.hat%*%replicate(T,as.vector(est_model_MCPanel_w$B)) + replicate(T,est_model_MCPanel_w$u) + t(replicate(N,est_model_MCPanel_w$v)) # use X with imputed endogenous values

    if(rev){
      est_model_MCPanel_w$impact <- (est_model_MCPanel_w$Mhat-Y)
      
      if(est_eastern){
        att <- as.matrix(colMeans(est_model_MCPanel_w$impact[,1:(t0-1)][rownames(est_model_MCPanel_w$impact) %in% eastern,])) # get mean pre-period impact on treated
        return(att)
      }
      if(est_swiss){
        att <- as.matrix(colMeans(est_model_MCPanel_w$impact[,1:(t0-1)][rownames(est_model_MCPanel_w$impact) %in% swiss,]))
        return(att)
      } else{
        return(est_model_MCPanel_w$impact)
      }
    }else{
      est_model_MCPanel_w$impact <- (Y-est_model_MCPanel_w$Mhat)
      if(est_eastern){
        att <- as.matrix(colMeans(est_model_MCPanel_w$impact[,t0:ncol(est_model_MCPanel_w$impact)][rownames(est_model_MCPanel_w$impact) %in% eastern,])) # get mean post-period impact on treated
        return(att)
      }
      if(est_swiss){
        as.matrix(colMeans(est_model_MCPanel_w$impact[,t0:ncol(est_model_MCPanel_w$impact)][rownames(est_model_MCPanel_w$impact) %in% swiss,]))
        return(att) 
      } else{
        return(est_model_MCPanel_w$impact)
      }
    }
  }else{
    ## ------
    ## MC-NNM
    ## ------
    
    est_model_MCPanel <- mcnnm(M = Y_obs, mask = treat_mat, W = weights, to_estimate_u = 1, to_estimate_v = 1, lambda_L=best_L, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
    
    est_model_MCPanel$Mhat <- est_model_MCPanel$L + replicate(T,est_model_MCPanel$u) + t(replicate(N,est_model_MCPanel$v))
    
    if(rev){
      est_model_MCPanel$impact <- (est_model_MCPanel$Mhat-Y)
      if(est_eastern){
        att <- as.matrix(colMeans(est_model_MCPanel$impact[,1:(t0-1)][rownames(est_model_MCPanel$impact) %in% eastern,])) # get mean pre-period impact on treated
        return(att)
      }
      if(est_swiss){
        att <- as.matrix(colMeans(est_model_MCPanel$impact[,1:(t0-1)][rownames(est_model_MCPanel$impact) %in% swiss,]))
        return(att)
      } else{
        return(est_model_MCPanel$impact)
      }
    }else{
      est_model_MCPanel$impact <- (Y-est_model_MCPanel$Mhat)
      if(est_eastern){
        att <- as.matrix(colMeans(est_model_MCPanel$impact[,t0:ncol(est_model_MCPanel$impact)][rownames(est_model_MCPanel$impact) %in% eastern,])) # get mean post-period impact on treated
        return(att)
      }
      if(est_swiss){
        att <- as.matrix(colMeans(est_model_MCPanel$impact[,t0:ncol(est_model_MCPanel$impact)][rownames(est_model_MCPanel$impact) %in% swiss,])) 
        return(att)
      } else{
        return(est_model_MCPanel$impact)
      }
    }
  }
}