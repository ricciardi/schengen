MCEst <- function(outcomes,rev=TRUE,covars=TRUE,fe=TRUE) {
  
  Y <- outcomes$M # NxT 
  
  treat <- outcomes$mask # NxT masked matrix 
  
  N <- nrow(treat)
  T <- ncol(treat)
  
  treat_mat <- 1-treat
  
  Y_obs <- Y * treat_mat
  
  W <- outcomes$W
  W <- W[rownames(W) %in% row.names(Y),]
  W <- W[row.names(Y),]  # ensure correct order
  
  z.cbw.eastern <- outcomes$z.cbw.eastern
  z.cbw.swiss <- outcomes$z.cbw.eastern
  
  weights <- matrix(NA, nrow=nrow(W), ncol=ncol(W), dimnames = list(rownames(W), colnames(W)))
  weights <- treat*(1-W) + (1-treat)*(W) 
  weights[rownames(weights) %in% outcomes$eastern,] <- weights[rownames(weights) %in% outcomes$eastern,] %*%diag(z.cbw.eastern)
  weights[rownames(weights) %in% outcomes$swiss,] <- weights[rownames(weights) %in% outcomes$swiss,] %*%diag(z.cbw.swiss)
  
  if(covars){
    
    X <- outcomes$X # NxT
    X.hat <- outcomes$X.hat # imputed endogenous values
    ## ------
    ## MC-NNM-W
    ## ------
    
    if(fe){
      est_model_MCPanel_w <- mcnnm_wc_cv(M = Y_obs, C = X, mask = treat_mat, W = weights, to_normalize = 1, to_estimate_u = 1, to_estimate_v = 1, num_lam_L = 10, num_lam_B = 5, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 2, is_quiet = 1) 
      est_model_MCPanel_w$Mhat <- est_model_MCPanel_w$L + X.hat%*%replicate(T,as.vector(est_model_MCPanel_w$B)) + replicate(T,est_model_MCPanel_w$u) + t(replicate(N,est_model_MCPanel_w$v)) # use X with imputed endogenous values
      est_model_MCPanel_w$rankL <- rankMatrix(t(est_model_MCPanel_w$L), method="qr.R")[1]
    }else{
      est_model_MCPanel_w <- mcnnm_wc_cv(M = Y_obs, C = X, mask = treat_mat, W = weights, to_normalize = 1, to_estimate_u = 0, to_estimate_v = 0, num_lam_L = 10, num_lam_B = 5, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 2, is_quiet = 1) 
      est_model_MCPanel_w$Mhat <- est_model_MCPanel_w$L + X.hat%*%replicate(T,as.vector(est_model_MCPanel_w$B))  # use X with imputed endogenous values
      est_model_MCPanel_w$rankL <- rankMatrix(t(est_model_MCPanel_w$L), method="qr.R")[1]
    }

    if(rev){
      est_model_MCPanel_w$impact <- (est_model_MCPanel_w$Mhat-Y)
    }else{
      est_model_MCPanel_w$impact <- (Y-est_model_MCPanel_w$Mhat)
    }
    
    return(est_model_MCPanel_w)
  } else{
    ## ------
    ## MC-NNM
    ## ------
    
    if(fe){
      est_model_MCPanel <- mcnnm_cv(M = Y_obs, mask = treat_mat, W = weights, to_estimate_u = 1, to_estimate_v = 1, num_lam_L = 10, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 2, is_quiet = 1)
      est_model_MCPanel$Mhat <- est_model_MCPanel$L + replicate(T,est_model_MCPanel$u) + t(replicate(N,est_model_MCPanel$v))
      est_model_MCPanel$rankL <- rankMatrix(t(est_model_MCPanel$L), method="qr.R")[1]
      }else{
        est_model_MCPanel <- mcnnm_cv(M = Y_obs, mask = treat_mat, W = weights, to_estimate_u = 0, to_estimate_v = 0, num_lam_L = 10, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 2, is_quiet = 1)
        est_model_MCPanel$Mhat <- est_model_MCPanel$L
        est_model_MCPanel$rankL <- rankMatrix(t(est_model_MCPanel$L), method="qr.R")[1]
      }

    if(rev){
      est_model_MCPanel$impact <- (est_model_MCPanel$Mhat-Y)
    }else{
      est_model_MCPanel$impact <- (Y-est_model_MCPanel$Mhat)
    }
    
    return(est_model_MCPanel)
  }
}