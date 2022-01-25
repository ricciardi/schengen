MCEst <- function(outcomes,cluster=c('eastern','swiss'),rev=TRUE,covars=TRUE,prop.model=FALSE,ADH=FALSE,DID=FALSE,IFE=FALSE) {
  
  if(cluster=='eastern'){
    Y <- outcomes$M[!rownames(outcomes$M)%in%outcomes$swiss,] # exclude Swiss regions
  }
  if(cluster=='swiss'){
    Y <- outcomes$M[!rownames(outcomes$M)%in%outcomes$eastern,] # exclude eastern regions
  }

  treat <- outcomes$mask
  treat <- treat[rownames(treat) %in% row.names(Y),]
  
  N <- nrow(treat)
  T <- ncol(treat)
  
  treat_mat <- 1-treat
  
  Y_obs <- Y * treat_mat
  
  W <- outcomes$W
  W <- W[rownames(W) %in% row.names(Y),]
  W <- W[row.names(Y),]  # ensure correct order
  
  ST <- intersect(rownames(outcomes$M)[outcomes$ST], row.names(Y))
  AT <- intersect(rownames(outcomes$M)[outcomes$AT], row.names(Y))
  
  z_weights <- outcomes$z_weights[rownames(outcomes$z_weights) %in% row.names(Y),]
  
  weights <- matrix(NA, nrow=nrow(W), ncol=ncol(W), dimnames = list(rownames(W), colnames(W)))

  weights[ST,] <- (1-diag(z_weights[ST,])*W[ST,])/(diag(z_weights[ST,])*W[ST,]) # elapsed-time weighting
  weights[AT,] <- (1-W[AT,])/(W[AT,]) # elapsed-time weighting
  
  if(covars){
    print("MC-NNM (weights + covars)")
    
    X <- outcomes$X # NxT
    X.hat <- outcomes$X.hat # imputed endogenous values
    
    X <- X[rownames(X) %in% row.names(Y),]
    X.hat <- X.hat[rownames(X.hat) %in% row.names(Y),]
    
    ## ------
    ## MC-NNM-W
    ## ------
    
    if(prop.model){
      est_model_MCPanel_w <- mcnnm_wc_cv(M = Y_obs, C = X.hat, mask = treat_mat, W = matrix(1, nrow(treat_mat),ncol(treat_mat)), to_normalize = 1, to_estimate_u = 1, to_estimate_v = 1, num_lam_L = 30, num_lam_B = 30, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 5, is_quiet = 1) # use X with imputed endogenous values # weights are equal
      est_model_MCPanel_w$Mhat <- plogis(est_model_MCPanel_w$L + X.hat%*%replicate(T,as.vector(est_model_MCPanel_w$B)) + replicate(T,est_model_MCPanel_w$u) + t(replicate(N,est_model_MCPanel_w$v))) # use X with imputed endogenous values
      est_model_MCPanel_w$rankL <- rankMatrix(t(est_model_MCPanel_w$L), method="qr.R")[1]
    }else{
      est_model_MCPanel_w <- mcnnm_wc_cv(M = Y_obs, C = X, mask = treat_mat, W = weights, to_normalize = 1, to_estimate_u = 1, to_estimate_v = 1, num_lam_L = 30, num_lam_B = 30, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 5, is_quiet = 1) 
      est_model_MCPanel_w$Mhat <- est_model_MCPanel_w$L + X.hat%*%replicate(T,as.vector(est_model_MCPanel_w$B)) + replicate(T,est_model_MCPanel_w$u) + t(replicate(N,est_model_MCPanel_w$v)) # use X with imputed endogenous values
      est_model_MCPanel_w$rankL <- rankMatrix(t(est_model_MCPanel_w$L), method="qr.R")[1]
    }

    if(rev){
      est_model_MCPanel_w$impact <- (est_model_MCPanel_w$Mhat-Y)
    }else{
      est_model_MCPanel_w$impact <- (Y-est_model_MCPanel_w$Mhat)
    }
    
    return(est_model_MCPanel_w)
  } else if(ADH){
    print("ADH")
    
    est_model_ADH <- list()
    est_model_ADH$Mhat <- adh_mp_rows(M=Y_obs, mask=treat_mat)
    
    if(rev){
      est_model_ADH$impact <- (est_model_ADH$Mhat-Y)
    }else{
      est_model_ADH$impact <- (Y-est_model_ADH$Mhat)
    }
    return(est_model_ADH)
  } else if(DID){
    print("DID")
    
    est_model_DID <- list()
    est_model_DID$Mhat <- DID(M=Y_obs, mask=treat_mat)
    
    if(rev){
      est_model_DID$impact <- (est_model_DID$Mhat-Y)
    }else{
      est_model_DID$impact <- (Y-est_model_DID$Mhat)
    }
    return(est_model_DID)
    
  } else if(IFE){
    print("IFE")
    
    est_model_IFE <- list()
    est_model_IFE$Mhat <- IFE(M=Y_obs, mask=treat_mat, k=2)
    
    if(rev){
      est_model_IFE$impact <- (est_model_IFE$Mhat-Y)
    }else{
      est_model_IFE$impact <- (Y-est_model_IFE$Mhat)
    }
    return(est_model_IFE)
    
  } else{
    print("MC-NNM (weights)")
    
    ## ------
    ## MC-NNM
    ## ------
    
    est_model_MCPanel <- mcnnm_cv(M = Y_obs, mask = treat_mat, W = weights, to_estimate_u = 1, to_estimate_v = 1, num_lam_L = 30, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 5, is_quiet = 1)
    est_model_MCPanel$Mhat <- est_model_MCPanel$L + replicate(T,est_model_MCPanel$u) + t(replicate(N,est_model_MCPanel$v))
    est_model_MCPanel$rankL <- rankMatrix(t(est_model_MCPanel$L), method="qr.R")[1]
    
    if(rev){
      est_model_MCPanel$impact <- (est_model_MCPanel$Mhat-Y)
    }else{
      est_model_MCPanel$impact <- (Y-est_model_MCPanel$Mhat)
    }
    
    return(est_model_MCPanel)
  }
}