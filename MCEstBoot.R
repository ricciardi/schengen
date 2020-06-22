MCEstBoot <- function(tseries,mask,W,treated,control,covars=TRUE) {
  
  Y <- t(tseries) # NxT 
  
  treat <- mask # NxT masked matrix 
  
  N <- nrow(treat)
  T <- ncol(treat)
  
  treat_mat <- 1-treat
  
  Y_obs <- Y * treat_mat
  
  if(covars){
    
    weights <- W
    weights <- weights[rownames(weights) %in% row.names(Y),]
    weights <- weights[row.names(Y),]  # reorder
    
    weights[rownames(weights) %in% control,] <- (weights[rownames(weights) %in% control,])/(1-weights[rownames(weights) %in% control,]) # control
    weights[rownames(weights) %in% treated,] <- 1/weights[rownames(weights) %in% treated,] # treated
    
    ## ------
    ## MC-NNM-W
    ## ------
    
    est_model_MCPanel_w <- mcnnm_wc_cv(M = Y_obs, C = weights, mask = treat_mat, W = weights, to_normalize = 1, to_estimate_u = 1, to_estimate_v = 1, num_lam_L = 10, num_lam_B = 10, niter = 1000, rel_tol = 1e-03, cv_ratio = 0.5, num_folds = 2, is_quiet = 1) 
    
    est_model_MCPanel_w$Mhat <- est_model_MCPanel_w$L + est_model_MCPanel_w$C*est_model_MCPanel_w$B + replicate(T,est_model_MCPanel_w$u) + t(replicate(N,est_model_MCPanel_w$v))
    
    est_model_MCPanel_w$impact <- (Y-est_model_MCPanel_w$Mhat)

    return(est_model_MCPanel_w$impact)
  } else{
    ## ------
    ## MC-NNM
    ## ------
    
    est_model_MCPanel <- mcnnm_cv(M = Y_obs, mask = treat_mat, W = weights, to_estimate_u = 1, to_estimate_v = 1,  num_lam_L = 10, niter = 1000, rel_tol = 1e-03, cv_ratio = 0.5, num_folds = 2, is_quiet = 1)
    
    est_model_MCPanel$Mhat <- est_model_MCPanel$L + replicate(T,est_model_MCPanel$u) + t(replicate(N,est_model_MCPanel$v))
    
    est_model_MCPanel$impact <- (Y-est_model_MCPanel$Mhat)
    
    return(est_model_MCPanel$impact)
  }
}