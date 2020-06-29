MCEst <- function(outcomes,rev=TRUE,covars=TRUE,nofes=FALSE) {
  
  Y <- outcomes$M # NxT 
  
  treat <- outcomes$mask # NxT masked matrix 
  
  N <- nrow(treat)
  T <- ncol(treat)
  
  treat_mat <- 1-treat
  
  Y_obs <- Y * treat_mat
  
  weights <- outcomes$W
  weights <- weights[rownames(weights) %in% row.names(Y),]
  weights <- weights[row.names(Y),]  # ensure correct order
  
  weights <- (weights)/(1-weights) 

  if(covars){
    
    X <- outcomes$X # NxT
    X.hat <- outcomes$X.hat # imputed endogenous values
    ## ------
    ## MC-NNM-W
    ## ------
    
    if(nofes){
      est_model_MCPanel_w <- mcnnm_wc_cv(M = Y_obs, C = X, mask = treat_mat, W = weights, to_normalize = 1, to_estimate_u = 0, to_estimate_v = 0, num_lam_L = 3, num_lam_B = 3, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 3, is_quiet = 1) 
      est_model_MCPanel_w$Mhat <- est_model_MCPanel_w$L + X.hat%*%replicate(T,as.vector(est_model_MCPanel_w$B)) # use X with imputed endogenous values
    }else{
      est_model_MCPanel_w <- mcnnm_wc_cv(M = Y_obs, C = X, mask = treat_mat, W = weights, to_normalize = 1, to_estimate_u = 1, to_estimate_v = 1, num_lam_L = 3, num_lam_B = 3, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 3, is_quiet = 1) 
      est_model_MCPanel_w$Mhat <- est_model_MCPanel_w$L + X.hat%*%replicate(T,as.vector(est_model_MCPanel_w$B)) + replicate(T,est_model_MCPanel_w$u) + t(replicate(N,est_model_MCPanel_w$v)) # use X with imputed endogenous values
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
    
    if(nofes){
      est_model_MCPanel <- mcnnm_cv(M = Y_obs, mask = treat_mat, W = weights, to_estimate_u = 0, to_estimate_v = 0, num_lam_L = 10, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 3, is_quiet = 1)
      est_model_MCPanel$Mhat <- est_model_MCPanel$L 
    }else{
      est_model_MCPanel <- mcnnm_cv(M = Y_obs, mask = treat_mat, W = weights, to_estimate_u = 1, to_estimate_v = 1, num_lam_L = 10, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 3, is_quiet = 1)
      est_model_MCPanel$Mhat <- est_model_MCPanel$L + replicate(T,est_model_MCPanel$u) + t(replicate(N,est_model_MCPanel$v))
    }

    if(rev){
      est_model_MCPanel$impact <- (est_model_MCPanel$Mhat-Y)
    }else{
      est_model_MCPanel$impact <- (Y-est_model_MCPanel$Mhat)
    }
    
    return(est_model_MCPanel)
  }
}