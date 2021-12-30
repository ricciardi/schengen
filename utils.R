##################################
# Utility fns.                  #
##################################

# mean exluding zero values (for calc. ATT)
nzmean <- function(x) {
  if (all(x==0)) 0 else mean(x[x!=0])
}

# convert NxT matrices to long panel

widetoLong <- function(M,mask,X=NULL,W=NULL){
  # inputs numeric matrices M, Mask, A, X, W
  # returns: data table with the columns
  ## person_id (numeric)
  ## year (character)
  ## y_value (numeric)
  ## W (numeric)
  ## x_value (numeric)
  ## w_value (numeric)
  
  data.long <- melt(M)
  colnames(data.long) <- c("person_id","year","y_value")
  
  data.long <- cbind(data.long, "W"=melt(mask)$value)
  
  if(!is.null(X)){
    data.long <- cbind(data.long, "x_value"=melt(X)$value)
  }
  if(!is.null(W)){
    data.long <- cbind(data.long, "w_value"=melt(W)$value)
  }
  
  data.long <-  data.long[order(data.long$person_id, data.long$year), ]
  
  return(data.table(data.long))
}

# revert long panel back to NxT matrices
longtoWide <- function(data.long){
  # inputs data table and outputs matrices M, mask, X W
  M.new <- as.matrix(reshape(data.long[,c("pid_boot","year_boot","y_value")], idvar = "pid_boot", timevar = "year_boot", direction = "wide")[,-1])
  mask.new <- as.matrix(reshape(data.long[,c("pid_boot","year_boot","W")], idvar = "pid_boot", timevar = "year_boot", direction = "wide")[,-1])
  if("x_value" %in% colnames(data.long)){
    X.new <- as.matrix(reshape(data.long[,c("pid_boot","year_boot","x_value")], idvar = "pid_boot", timevar = "year_boot", direction = "wide")[,-1])
  }else{
    X.new <- NULL
  }
  if("w_value" %in% colnames(data.long)){
    W.new <- as.matrix(reshape(data.long[,c("pid_boot","year_boot","w_value")], idvar = "pid_boot", timevar = "year_boot", direction = "wide")[,-1])
  }else{
    W.new <- NULL
  }
  return(list("M"=M.new, "mask"=mask.new, "X"=X.new, "W"=W.new))
}

#  one bootstrap sample
one_boot <- function(sim_num, current_data_realized_long,Y, estimator, z_weights){
  boot_matrices <- list()
  boot_matrices$mask <- matrix(0, nrow=nrow(Y), ncol=ncol(Y)) # treat matrix
  while(any(rowSums(boot_matrices$mask)<1) || max(rowSums(boot_matrices$mask))<T){ # ensure that there are LT and AT units
    num_units <- data.table::uniqueN(current_data_realized_long$person_id)
    sample_units <- data.table::data.table((table(sample(1:num_units, replace =  TRUE))))
    sample_units[, person_id := as.numeric(V1)]
    sample_units[, N := as.numeric(N)]
    
    boot_DT <- merge(current_data_realized_long, sample_units,
                     by = "person_id", all.x = TRUE)
    boot_DT <- boot_DT[!is.na(N),]
    boot_DT[, ID := .I]
    
    boot_DT <- boot_DT[rep(boot_DT$ID, boot_DT$N)]
    boot_DT$pid_boot <- rep(1:nrow(Y), each = nrow(Y))
    boot_DT$year_boot <- rep(1:ncol(Y), times = ncol(Y))
    
    boot_matrices <- longtoWide(data.long=boot_DT) 
    
    ST.boot <- which(rowSums(boot_matrices$mask)<ncol(Y)) # switch treated indices
  }
  
  if(!is.null(z_weights)){
    # Impute endogenous covariate values
    
    boot_model_endogenous <- mcnnm(M = boot_matrices$X, mask = boot_matrices$mask, W = matrix(1, nrow(boot_matrices$mask),ncol(boot_matrices$mask)), to_estimate_u = 1, to_estimate_v = 1, 
                                   lambda_L = 0.02146935, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]] # outcome is covariate, prop. weights are equal
    boot_model_endogenous$Mhat <- boot_model_endogenous$L + replicate(ncol(boot_matrices$mask),boot_model_endogenous$u) + t(replicate(nrow(boot_matrices$mask),boot_model_endogenous$v))
    
    X.hat.boot <- boot_matrices$X*boot_matrices$mask + boot_model_endogenous$Mhat*(1-boot_matrices$mask) # only endogenous values imputed
    
    # Estimate propensity weights by matrix completion
    
    boot_model_pweights <- mcnnm_wc(M = (1-boot_matrices$mask), C = X.hat.boot, mask = matrix(1, nrow(boot_matrices$mask),ncol(boot_matrices$mask)), # no missing entries
                                    W = matrix(1, nrow(boot_matrices$mask),ncol(boot_matrices$mask)), to_normalize = 1, to_estimate_u = 1, to_estimate_v = 1, lambda_L = 9.29274e-05, lambda_B = 0.00886265, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]] # use X with imputed endogenous values
    boot_model_pweights$Mhat <- plogis(boot_model_pweights$L + X.hat.boot%*%replicate(ncol(boot_matrices$mask),as.vector(boot_model_pweights$B)) + replicate(ncol(boot_matrices$mask),boot_model_pweights$u) + t(replicate(nrow(boot_matrices$mask), boot_model_pweights$v)))
    
    weights.boot <- (1-diag(z_weights)*boot_model_pweights$Mhat)/(diag(z_weights)*boot_model_pweights$Mhat) # elapsed-time weighting
  }
  
  # estimators
  if(estimator=="mc_plain"){
    boot_mc_plain <- mcnnm(M = boot_matrices$M, mask = boot_matrices$mask, W = matrix(1, nrow(boot_matrices$mask),ncol(boot_matrices$mask)), to_estimate_u = 1, to_estimate_v = 1, lambda_L = 0.0014828, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
    boot_mc_plain$Mhat <- boot_mc_plain$L + replicate(ncol(boot_matrices$M),boot_mc_plain$u) + t(replicate(nrow(boot_matrices$M),boot_mc_plain$v))
    tau.boot <- (boot_mc_plain$Mhat-Y) # estimated treatment effect
  }
  
  if(estimator=="mc_covars"){
    boot_mc_covars <- mcnnm_wc(M = boot_matrices$M, C = boot_matrices$X, mask =  boot_matrices$mask, W = matrix(1, nrow( boot_matrices$mask),ncol(boot_matrices$mask)), to_estimate_u = 1, to_estimate_v = 1, lambda_L = 0.00108901, lambda_B = 0.100334, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
    boot_mc_covars$Mhat <- boot_mc_covars$L + X.hat.boot%*%replicate(ncol(boot_matrices$mask),as.vector(boot_mc_covars$B)) + replicate(ncol(boot_matrices$mask),boot_mc_covars$u) + t(replicate(nrow(boot_matrices$mask),boot_mc_covars$v))
    tau.boot <- (boot_mc_covars$Mhat-Y) # estimated treatment effect
  }
  if(estimator=="mc_weights"){
    boot_mc_weights <- mcnnm(M = boot_matrices$M, mask = boot_matrices$mask, W = weights.boot, to_estimate_u = 1, to_estimate_v = 1, lambda_L = 0.00333718, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
    boot_mc_weights$Mhat <- boot_mc_weights$L + replicate(ncol(boot_matrices$mask),boot_mc_weights$u) + t(replicate(nrow(boot_matrices$mask),boot_mc_weights$v))
    tau.boot <- (boot_mc_weights$Mhat-Y) # estimated treatment effect
  }
  if(estimator=="mc_weights_covars"){
    boot_mc_weights_covars <- mcnnm_wc(M = boot_matrices$M, C = boot_matrices$X, mask = boot_matrices$mask, W = weights.boot, to_estimate_u = 1, to_estimate_v = 1, lambda_L = 0.00108901, lambda_B = 0.100334, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
    boot_mc_weights_covars$Mhat <- boot_mc_weights_covars$L + X.hat.boot%*%replicate(ncol(boot_matrices$mask),as.vector(boot_mc_weights_covars$B)) + replicate(ncol(boot_matrices$mask),boot_mc_weights_covars$u) + t(replicate(nrow(boot_matrices$mask),boot_mc_weights_covars$v))
    tau.boot <- (boot_mc_weights_covars$Mhat-Y) # estimated treatment effect
  }
  if(estimator=="ADH"){
    boot_model_ADH <- list()
    boot_model_ADH$Mhat <- adh_mp_rows(M=boot_matrices$M, mask=boot_matrices$mask)
    tau.boot <- (boot_model_ADH$Mhat-Y) # estimated treatment effect
  }
  if(estimator=="DID"){
    boot_model_DID <- list()
    boot_model_DID$Mhat <- DID(M=boot_matrices$M, mask=boot_matrices$mask)
    tau.boot <- (boot_model_DID$Mhat-Y) # estimated treatment effect
  }
  
  if(estimator=="IFE"){
    boot_model_IFE <- list()
    boot_model_IFE$Mhat <- IFE(M=boot_matrices$M, mask=boot_matrices$mask)
    tau.boot <- (boot_model_IFE$Mhat-Y) # estimated treatment effect
  }
  
  # Calc. real ATT on the ST
  att.boot <- mean(apply(tau.boot*(1-boot_matrices$mask),1,nzmean)[ST.boot])
  return(att.boot)
}

# bootstap with cluster
clustered_bootstrap <- function(current_data_realized_long, estimator, Y, B = 999, z_weights=NULL){
  clustered_bootstrap_var <- var(unlist(lapply(c(1:B), one_boot, current_data_realized_long, Y, estimator, z_weights)))
  return(clustered_bootstrap_var)
}

# confidence interval
CI_test <- function(est_coefficent, real_coefficent, est_var,alpha=0.05){
  as.numeric(est_coefficent - qnorm(1 - alpha/2)*sqrt(est_var) <= real_coefficent &
               est_coefficent + qnorm(1 - alpha/2)*sqrt(est_var) >= real_coefficent )
}