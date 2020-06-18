MCEstBoot <- function(tseries,mask,W,covars=NULL) {
  
  Y <- t(tseries) # NxT 

  weights <- W
  weights <- weights[rownames(weights) %in% row.names(Y),]
  weights <- weights[row.names(Y),]  # reorder
  
  treat <- mask # NxT masked matrix 
  
  N <- nrow(treat)
  T <- ncol(treat)
  
  treat_mat <- 1-treat
  
  Y_obs <- Y * treat_mat
  
  if(!is.null(covars)){
    ## ------
    ## MC-NNM-W
    ## ------
    
    lam.min <- mcnnm_wc_lam_range(M=Y_obs, C=(weights)/(1-weights), X=matrix(0L,0,0), Z=matrix(0L,0,0), mask=treat_mat, W=(weights)/(1-weights), 
                                  to_estimate_u = 1L, to_estimate_v = 1L, niter = 1000L, rel_tol = 1e-05)
    
    est_model_MCPanel_w <- mcnnm_wc_fit(M=Y_obs, C=(weights)/(1-weights), X=matrix(0L,0,0), Z=matrix(0L,0,0), mask=treat_mat, W=(weights)/(1-weights),
                                        lambda_L=lam.min$lambda_L_max, lambda_H=lam.min$lambda_H_max, lambda_B=lam.min$lambda_B_max, 
                                        to_normalize = 1, to_estimate_u = 1, to_estimate_v = 1, to_add_ID = 1, niter = 1000, rel_tol = 1e-05, is_quiet = 1) 
    
    est_model_MCPanel_w$Mhat <- est_model_MCPanel_w$L + est_model_MCPanel_w$C%*%est_model_MCPanel_w$B + replicate(T,est_model_MCPanel_w$u) + t(replicate(N,est_model_MCPanel_w$v))
    
    est_model_MCPanel_w$impact <- (Y-est_model_MCPanel_w$Mhat)

    return(est_model_MCPanel_w$impact)
  } else{
    ## ------
    ## MC-NNM
    ## ------
    
    lam.min <- mcnnm_lam_range(Y_obs, treat_mat, (weights)/(1-weights), to_estimate_u = 1L, to_estimate_v = 1L,
                               niter = 1000L, rel_tol = 1e-05)
    
    est_model_MCPanel <- mcnnm_fit(Y_obs, treat_mat, (weights)/(1-weights), to_estimate_u = 1, to_estimate_v = 1, 
                                   lambda_L = lam.min, niter = 1000, rel_tol = 1e-05, is_quiet = 1)
    
    est_model_MCPanel$Mhat <- est_model_MCPanel$L + replicate(T,est_model_MCPanel$u) + t(replicate(N,est_model_MCPanel$v))

    est_model_MCPanel$impact <- (Y-est_model_MCPanel$Mhat)
    
    return(est_model_MCPanel$impact)
  }
}