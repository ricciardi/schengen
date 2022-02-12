##################################
# Utility fns.                  #
##################################

# function to bound probabilities to be used when making predictions
boundProbs <- function(x,bounds=c(0,1)){
  x[x>max(bounds)] <- max(bounds)
  x[x<min(bounds)] <- min(bounds)
  return(x)
}

DropVariance <- function(mat){ # Remove units with no variance
  drop <- names(which(apply(t(mat), 2, var) == 0))
  return(as.matrix(mat[!rownames(mat)%in%drop,]))
}

# mean exluding zero values (for calc. ATT)
nzmean <- function(x) {
  if (all(x==0)) 0 else mean(x[x!=0])
}

# convert NxT matrices to long panel

widetoLong <- function(Y,mask,X=NULL){
  # inputs numeric matrices Y, Mask, X
  # returns: data table with the columns

  
  data.long <- melt(Y)
  colnames(data.long) <- c("person_id","year","y_value")
  
  data.long <- cbind(data.long, "W"=melt(mask)$value)
  
  if(!is.null(X)){
    data.long <- cbind(data.long, "x_value"=melt(X)$value)
  }

  data.long <-  data.long[order(data.long$person_id, data.long$year), ]
  
  return(data.table(data.long))
}

# revert long panel back to NxT matrices
longtoWide <- function(data.long){
  # inputs data table and outputs matrices Y, mask, X
  Y.new <- as.matrix(reshape(data.long[,c("pid_boot","year_boot","y_value")], idvar = "pid_boot", timevar = "year_boot", direction = "wide")[,-1])
  mask.new <- as.matrix(reshape(data.long[,c("pid_boot","year_boot","W")], idvar = "pid_boot", timevar = "year_boot", direction = "wide")[,-1])
  if("x_value" %in% colnames(data.long)){
    X.new <- as.matrix(reshape(data.long[,c("pid_boot","year_boot","x_value")], idvar = "pid_boot", timevar = "year_boot", direction = "wide")[,-1])
  }else{
    X.new <- NULL
  }
  return(list("Y"=Y.new, "mask"=mask.new, "X"=X.new))
}

#  one bootstrap sample
one_boot <- function(sim_num, current_data_realized_long, N, T,estimator, est_weights, return_per_period=FALSE, placebo=FALSE){
  boot_matrices <- list()
  boot_matrices$mask <- matrix(0, nrow=N, ncol=T) # treat matrix
  while(any(rowSums(boot_matrices$mask)<=1) || max(rowSums(boot_matrices$mask))<T){ # ensure that there are ST (switch after at least 1 period) and AT units
    num_units <- data.table::uniqueN(current_data_realized_long$person_id)
    sample_units <- data.table::data.table((table(sample(1:num_units, replace =  TRUE))))
    sample_units[, person_id := as.numeric(V1)]
    sample_units[, N := as.numeric(N)]
    
    boot_DT <- merge(current_data_realized_long, sample_units,
                     by = "person_id", all.x = TRUE)
    boot_DT <- boot_DT[!is.na(N),]
    boot_DT[, ID := .I]
    
    boot_DT <- boot_DT[rep(boot_DT$ID, boot_DT$N)]
    boot_DT$pid_boot <- rep(1:N, each = T)
    boot_DT$year_boot <- rep(1:T, times = N)
    
    boot_matrices <- longtoWide(data.long=boot_DT) 
  }
  
  # get vector of initial treatment periods
  
  A.boot <- aggregate(col ~ row,
                 data = which(boot_matrices$mask == 1, arr.ind = T),
                 FUN = function(x) x[1])$col
  A.boot[which(A.boot==1)] <- Inf
  
  ST.boot <- which(!is.infinite(A.boot)) # which(rowSums(boot_matrices$mask)<T) # switch treated indices 
  AT.boot <- which(is.infinite(A.boot)) # always-treated indices 
  
  # Observed outcome matrix 
  
  boot_matrices$obs_mat <- boot_matrices$Y * boot_matrices$mask
  
  if(est_weights){
    
    # Estimate z_weights
    boot_matrices$z_weights <- matrix(0,N,T,byrow = TRUE)

    for(i in ST.boot){
      boot_matrices$z_weights[i,] <- c(plogis(A.boot[i]:1, scale=8),plogis(1:(ncol(boot_matrices$mask)-A.boot[i]),scale=8))
    }
      
    # Impute endogenous covariate values
    
    boot_model_endogenous <- mcnnm(M = boot_matrices$X, mask = boot_matrices$mask, W = matrix(1, nrow(boot_matrices$mask),ncol(boot_matrices$mask)), to_estimate_u = 1, to_estimate_v = 1, 
                                   lambda_L = 0.02146935, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]] # outcome is covariate, prop. weights are equal
    boot_model_endogenous$Mhat <- boot_model_endogenous$L + replicate(ncol(boot_matrices$mask),boot_model_endogenous$u) + t(replicate(nrow(boot_matrices$mask),boot_model_endogenous$v))
    
    X.hat.boot <- boot_matrices$X*boot_matrices$mask + boot_model_endogenous$Mhat*(1-boot_matrices$mask) # only endogenous values imputed
    
    # Estimate propensity weights by matrix completion
    
    boot_model_pweights <- mcnnm_wc(M = (1-boot_matrices$mask), C = X.hat.boot, mask = matrix(1, nrow(boot_matrices$mask),ncol(boot_matrices$mask)), # no missing entries
                                    W = matrix(1, nrow(boot_matrices$mask),ncol(boot_matrices$mask)), to_normalize = 1, to_estimate_u = 1, to_estimate_v = 1, lambda_L = 9.29274e-05, lambda_B = 0.00886265, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]] # use X with imputed endogenous values
    boot_model_pweights$Mhat <- boot_model_pweights$L + X.hat.boot%*%replicate(ncol(boot_matrices$mask),as.vector(boot_model_pweights$B)) + replicate(ncol(boot_matrices$mask),boot_model_pweights$u) + t(replicate(nrow(boot_matrices$mask), boot_model_pweights$v))
    boot_model_pweights$Mhat <- boundProbs(boot_model_pweights$Mhat)
      
    weights.boot <-  matrix(0, nrow=N, ncol=T) # treat matrix
    
    weights.boot[ST.boot,] <- (1-diag(boot_matrices$z_weights[ST.boot,])*boot_model_pweights$Mhat[ST.boot,])/(diag(boot_matrices$z_weights[ST.boot,])*boot_model_pweights$Mhat[ST.boot,]) # elapsed-time weighting
    weights.boot[AT.boot,] <- (1-boot_model_pweights$Mhat[AT.boot,])/(boot_model_pweights$Mhat[AT.boot,]) # elapsed-time weighting
  }
  
  # estimators
  if(estimator=="mc_plain"){
    boot_mc_plain <- mcnnm(M = boot_matrices$obs_mat, mask = boot_matrices$mask, W = matrix(1, nrow(boot_matrices$mask),ncol(boot_matrices$mask)), to_estimate_u = 1, to_estimate_v = 1, lambda_L = 0.0014828, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
    boot_mc_plain$Mhat <- boot_mc_plain$L + replicate(ncol(boot_matrices$obs_mat),boot_mc_plain$u) + t(replicate(nrow(boot_matrices$obs_mat),boot_mc_plain$v))
    tau.boot <- (boot_mc_plain$Mhat-boot_matrices$Y) # estimated treatment effect
  }
  
  if(estimator=="mc_covars"){
    boot_mc_covars <- mcnnm_wc(M = boot_matrices$obs_mat, C = boot_matrices$X, mask =  boot_matrices$mask, W = matrix(1, nrow( boot_matrices$mask),ncol(boot_matrices$mask)), to_estimate_u = 1, to_estimate_v = 1, lambda_L = 0.00108901, lambda_B = 0.100334, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
    boot_mc_covars$Mhat <- boot_mc_covars$L + X.hat.boot%*%replicate(ncol(boot_matrices$mask),as.vector(boot_mc_covars$B)) + replicate(ncol(boot_matrices$mask),boot_mc_covars$u) + t(replicate(nrow(boot_matrices$mask),boot_mc_covars$v))
    tau.boot <- (boot_mc_covars$Mhat-boot_matrices$Y) # estimated treatment effect
  }
  if(estimator=="mc_weights"){
    boot_mc_weights <- mcnnm(M = boot_matrices$obs_mat, mask = boot_matrices$mask, W = weights.boot, to_estimate_u = 1, to_estimate_v = 1, lambda_L = 0.00333718, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
    boot_mc_weights$Mhat <- boot_mc_weights$L + replicate(ncol(boot_matrices$mask),boot_mc_weights$u) + t(replicate(nrow(boot_matrices$mask),boot_mc_weights$v))
    tau.boot <- (boot_mc_weights$Mhat-boot_matrices$Y) # estimated treatment effect
  }
  if(estimator=="mc_weights_covars"){
    boot_mc_weights_covars <- mcnnm_wc(M = boot_matrices$obs_mat, C = boot_matrices$X, mask = boot_matrices$mask, W = weights.boot, to_estimate_u = 1, to_estimate_v = 1, lambda_L = 0.00108901, lambda_B = 0.100334, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
    while(any(is.na(boot_mc_weights_covars$B))){
      boot_mc_weights_covars <- mcnnm_wc(M = boot_matrices$obs_mat, C = boot_matrices$X, mask = boot_matrices$mask, W = weights.boot, to_estimate_u = 1, to_estimate_v = 1, lambda_L = 0.00108901, lambda_B = 0.100334, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
    }
    boot_mc_weights_covars$Mhat <- boot_mc_weights_covars$L + X.hat.boot%*%replicate(ncol(boot_matrices$mask),as.vector(boot_mc_weights_covars$B)) + replicate(ncol(boot_matrices$mask),boot_mc_weights_covars$u) + t(replicate(nrow(boot_matrices$mask),boot_mc_weights_covars$v))
    tau.boot <- (boot_mc_weights_covars$Mhat-boot_matrices$Y) # estimated treatment effect
  }
  if(estimator=="ADH"){
    boot_model_ADH <- list()
    boot_model_ADH$Mhat <- adh_mp_rows(M=boot_matrices$obs_mat, mask=boot_matrices$mask)
    tau.boot <- (boot_model_ADH$Mhat-boot_matrices$Y) # estimated treatment effect
  }
  if(estimator=="DID"){
    boot_model_DID <- list()
    boot_model_DID$Mhat <- DID(M=boot_matrices$obs_mat, mask=boot_matrices$mask)
    tau.boot <- (boot_model_DID$Mhat-boot_matrices$Y) # estimated treatment effect
  }
  
  if(estimator=="IFE"){
    boot_model_IFE <- list()
    boot_model_IFE$Mhat <- IFE(M=boot_matrices$obs_mat, mask=boot_matrices$mask)
    tau.boot <- (boot_model_IFE$Mhat-boot_matrices$Y) # estimated treatment effect
  }
  
  # Calc. real ATT on the ST
  att.boot <- apply(tau.boot*(1-boot_matrices$mask),1,nzmean)[ST.boot]
  att.bar.boot <- mean(att.boot)
  if(return_per_period){
    return(list("att.bar.boot"=att.bar.boot, "att.boot"=att.boot))
  }else{
    return(att.bar.boot)
  }
}

# bootstap with cluster (return att.bar variance)
clustered_bootstrap <- function(current_data_realized_long, N,T,estimator, B = 999, est_weights, return_per_period=FALSE, return_replicates=FALSE, ncores=NULL){
  if(return_per_period){
    if(!is.null(ncores)){
      boot_stats <- mclapply(c(1:B), one_boot, current_data_realized_long, N, T, estimator, est_weights, return_per_period=TRUE, mc.cores=ncores)
    }else{
      boot_stats <- lapply(c(1:B), one_boot, current_data_realized_long, N, T, estimator, est_weights, return_per_period=TRUE)
    }
    clustered_bootstrap_att_var <- colVars(as.matrix(plyr::ldply(sapply(boot_stats, function(x) rbind(x$att.boot)), rbind)), na.rm = TRUE)
    clustered_bootstrap_att_bar_var <- var(unlist(sapply(boot_stats, function(x) c(x$att.bar.boot))))
    return(list("att.bar.boot.var"=clustered_bootstrap_att_bar_var,"att.boot.var"=clustered_bootstrap_att_var))
  }
  if(return_replicates){
    boot_stats <- unlist(lapply(c(1:B), one_boot, current_data_realized_long, N, T, estimator, est_weights, return_per_period=FALSE))
    return(boot_stats)
  }else{
    clustered_bootstrap_var <- var(unlist(lapply(c(1:B), one_boot, current_data_realized_long, N, T, estimator, est_weights, return_per_period=FALSE)))
    return(clustered_bootstrap_var)
  }
}

# confidence interval
CI_test <- function(est_coefficent, real_coefficent, est_var,alpha=0.05){
  return(as.numeric(est_coefficent - qnorm(1 - alpha/2)*sqrt(est_var) <= real_coefficent &
               est_coefficent + qnorm(1 - alpha/2)*sqrt(est_var) >= real_coefficent ))
}

boot_CI <- function(est_coefficent,est_var,alpha=0.05){
  return(list("lb"=est_coefficent - qnorm(1 - alpha/2)*sqrt(est_var),
               "ub"=est_coefficent + qnorm(1 - alpha/2)*sqrt(est_var) ))
}