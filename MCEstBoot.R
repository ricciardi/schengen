MCEstBoot <- function(tseries,mask,W=NULL,treated,control,covars=FALSE) {
  
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
    
    est_model_MCPanel_w <- mcnnm_wc_fit(M = Y_obs, C = weights, mask = treat_mat, W = weights, lambda_L=0.1, lambda_B=0.1, to_normalize = 1, to_estimate_u = 1, to_estimate_v = 1,
                                        niter = 1000, rel_tol = 1e-05, is_quiet = 1) 
    
    est_model_MCPanel_w$Mhat <- est_model_MCPanel_w$L + est_model_MCPanel_w$C*est_model_MCPanel_w$B + replicate(T,est_model_MCPanel_w$u) + t(replicate(N,est_model_MCPanel_w$v))
    
    est_model_MCPanel_w$impact <- (Y-est_model_MCPanel_w$Mhat)

    return(est_model_MCPanel_w$impact)
  } else{
    ## ------
    ## MC-NNM
    ## ------
    
    est_model_MCPanel <- mcnnm_fit(M = Y_obs, mask = treat_mat, lambda_L = 0.1, to_estimate_u = 1, to_estimate_v = 1, 
                                   niter = 1000, rel_tol = 1e-05, is_quiet = 1)
    
    est_model_MCPanel$Mhat <- est_model_MCPanel$L + replicate(T,est_model_MCPanel$u) + t(replicate(N,est_model_MCPanel$v))
    
    est_model_MCPanel$impact <- (Y-est_model_MCPanel$Mhat)
    
    return(est_model_MCPanel$impact)
  }
}