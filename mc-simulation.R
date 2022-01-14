##################################
# Matrix completion simulations  #
##################################

library(MCPanel)
library(matrixStats)
library(Matrix)
library(tictoc)
library(MASS)
library(data.table)
library(reshape)
library(reshape2)
library(emfactor)

source('utils.R')
source('IFE.R')

# Setup parallel processing
doMPI <- TRUE
if(doMPI){
  library(doMPI)
  
  # Start cluster
  cl <- startMPIcluster()
  
  # Register cluster
  registerDoMPI(cl)
  
  # Check cluster size
  print(paste0("cluster size: ", clusterSize(cl)))
  
} else{
  library(parallel)
  library(doParallel)
  library(foreach)
  
  cores <- parallel::detectCores()
  print(paste0("number of cores used: ", cores))
  
  cl <- parallel::makeCluster(cores, outfile="")
  
  doParallel::registerDoParallel(cl) # register cluster
}

MCsim <- function(N,T,R,T0,N_t,noise_sc,delta_sc,gamma_sc,beta_sc,shift_sc,n){
  
  # check inputs 
  if(N!=T){
    stop("Matrix must be square, N=T")
  }
  
  # set the seed
  print(paste0("run number: ", n))
  set.seed(n, "L'Ecuyer-CMRG")
  
  mask <- matrix(0, nrow=N, ncol=T) # treat matrix
  
  while(any(rowSums(mask)<1) || max(rowSums(mask))<T){ # ensure that there are LT and AT units
    
    # make the mean and var of the potential outcomes matrix
    
    mean_vec_0_N <- rep(0,N)
    mean_vec_1_N <- rep(1,N)
    
    sigma_mat_N <- diag(N)
    sigma_mat_N[sigma_mat_N==0] <- 0.2
    
    mean_vec_1_R <- rep(1, R)
    sigma_mat_R <- diag(R)
    sigma_mat_R[sigma_mat_R==0] <- 0.2
    
    # Create Matrices
    A <- mvrnorm(T, mu=mean_vec_1_R, Sigma=sigma_mat_R)  #replicate(T,rnorm(R))
    B <- mvrnorm(R, mu=mean_vec_1_N, Sigma=sigma_mat_N)  # replicate(R,rnorm(N))
    X <- mvrnorm(T, mu=mean_vec_1_N, Sigma=sigma_mat_N) # replicate(T,rnorm(N))
    delta <- delta_sc*mvrnorm(mu=mean_vec_0_N, Sigma=sigma_mat_N) #delta_sc*rnorm(N)
    gamma <- gamma_sc*mvrnorm(mu=mean_vec_0_N, Sigma=sigma_mat_N) #gamma_sc*rnorm(T)
    beta <- beta_sc*mvrnorm(mu=mean_vec_0_N, Sigma=sigma_mat_N) #beta_sc*rnorm(T)
    noise <- noise_sc*mvrnorm(T, mean_vec_0_N, Sigma=sigma_mat_N) # noise_sc*replicate(T,rnorm(N))
    
    # True outcome model
    true_mat <- A%*%B + X%*%replicate(T,as.vector(beta)) + replicate(T,delta) + t(replicate(N,gamma)) # potential outcomes under AT
    
    # True treatment model
    
    mean_vec_0_R <- rep(0,R)
    
    sigma_mat_N_treat <- diag(N)
    sigma_mat_N_treat[sigma_mat_N_treat==0] <- 0.3
    
    sigma_mat_R_treat <- diag(R)
    sigma_mat_R[sigma_mat_R_treat==0] <- 0.3
    
    C <-  mvrnorm(T, mu=mean_vec_0_R, Sigma=sigma_mat_R_treat)  #replicate(T,rnorm(R))
    D <- mvrnorm(R, mu=mean_vec_0_N, Sigma=sigma_mat_N_treat)  # replicate(R,rnorm(N))
    xi <- delta_sc*mvrnorm(mu=mean_vec_0_N, Sigma=sigma_mat_N) 
    psi <- gamma_sc*mvrnorm(mu=mean_vec_0_N, Sigma=sigma_mat_N) 
    phi <- beta_sc*mvrnorm(mu=mean_vec_0_N, Sigma=sigma_mat_N) 
    
    e <-plogis(C%*%D + X%*%replicate(T,as.vector(phi)) + replicate(T,xi) + t(replicate(N,psi))) # prob of being 0 (ST in pre-treatment)
    treat_mat <- stag_adapt(M = matrix(1, nrow=N, ncol=T) , N_t = N_t, T0 = T0, treat_indices = 0, weights = e[,T0]) # 0s missing and to be imputed (LT); 1s are observed (AT)
    
    mask <- treat_mat[,c(T:1)] # retrospective analysis
  }
  
  fr_obs <- sum(mask)/(N*T) # store fraction observed entries
  print(paste0("fraction observed: ", fr_obs))
  
  # get vector of initial treatment periods
  
  A <- aggregate(col ~ row,
                  data = which(mask == 1, arr.ind = T),
                  FUN = function(x) x[1])$col
  A[which(A==1)] <- Inf
  
  ST <- which(!is.infinite(A)) # switch treated indices
  AT <- which(is.infinite(A)) # always-treated indices
  
  # Elapsed time weighted treatment
  
  z_weights <- matrix(0,N,T,byrow = TRUE)
  
  for(i in ST){
    z_weights[i,] <- c(plogis(A[i]:1, scale=8),plogis(1:(ncol(mask)-A[i]),scale=8))
  }
  
  shift <- z_weights*shift_sc*(1-mask) 
  
  shifted_mat <- true_mat - shift # Y(AT) - Y(ST)
  
  # Calc. real ATT on the ST
  
  att.true <- mean(apply(shift,1,nzmean)[ST]) # the avg. of ATTs (unit-level)
  
  # Observed outcome matrix (+ noise)
  
  noisy_mat <- shifted_mat + noise # add noise
  
  obs_mat <- noisy_mat * mask 
  
  # Impute endogenous covariate values
  
  est_model_endogenous <- mcnnm(M = X, mask = mask, W = matrix(1, nrow(mask),ncol(mask)), to_estimate_u = 1, to_estimate_v = 1, 
                                lambda_L = 0.02146935, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]] # outcome is covariate, prop. weights are equal
  est_model_endogenous$Mhat <- est_model_endogenous$L + replicate(T,est_model_endogenous$u) + t(replicate(N,est_model_endogenous$v))
  
  X.hat <- X*mask + est_model_endogenous$Mhat*(1-mask) # only endogenous values imputed
  
  # Estimate propensity weights by matrix completion
  
  est_model_pweights <- mcnnm_wc(M = (1-mask), C = X.hat, mask =  matrix(1, nrow(mask),ncol(mask)), # no missing entries
                                 W = matrix(1, nrow(mask),ncol(mask)), to_normalize = 1, to_estimate_u = 1, to_estimate_v = 1, lambda_L = 9.29274e-05, lambda_B = 0.00886265, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]] # use X with imputed endogenous values
  est_model_pweights$Mhat <- plogis(est_model_pweights$L + X.hat%*%replicate(T,as.vector(est_model_pweights$B)) + replicate(T,est_model_pweights$u) + t(replicate(N,est_model_pweights$v)))
  
  weights <-  matrix(0, nrow=N, ncol=T) # treat matrix
    
  weights[ST,] <- (1-diag(z_weights[ST,])*est_model_pweights$Mhat[ST,])/(diag(z_weights[ST,])*est_model_pweights$Mhat[ST,]) # elapsed-time weighting
  weights[AT,] <- (1-est_model_pweights$Mhat[AT,])/(est_model_pweights$Mhat[AT,]) # elapsed-time weighting
  
  ## ------ ------ ------ ------ ------
  ## MC-NNM plain (no weighting, no covariate)
  ## ------ ------ ------ ------ ------
  
  est_mc_plain <- mcnnm(M = obs_mat, mask = mask, W = matrix(1, nrow(mask),ncol(mask)), to_estimate_u = 1, to_estimate_v = 1, lambda_L = 0.0014828, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
  est_mc_plain$Mhat <- est_mc_plain$L + replicate(T,est_mc_plain$u) + t(replicate(N,est_mc_plain$v))
  est_mc_plain$tau <- (est_mc_plain$Mhat-noisy_mat) # estimated treatment effect, Y(AT) - Y(ST)
  est_mc_plain$err <- (est_mc_plain$tau  - shift) # error (wrt to ground truth)
  
  est_mc_plain$msk_err <- est_mc_plain$err*(1-mask) # masked error (wrt to ground truth)
  est_mc_plain$RMSE <- sqrt((1/sum(1-mask)) * sum(est_mc_plain$msk_err^2)) # RMSE on test set (wrt to ground truth)
  print(paste("MC-NNM (Plain) RMSE:", round(est_mc_plain$RMSE,3)))
  
  est_mc_plain$att <- apply(est_mc_plain$tau*(1-mask),1,nzmean)[ST]
  est_mc_plain$att.bar <- mean(est_mc_plain$att)
  est_mc_plain$abs.bias <- abs(est_mc_plain$att.bar-att.true)
  est_mc_plain$rel.abs.bias <- est_mc_plain$abs.bias/att.true
  print(paste("MC-NNM (Plain) rel. abs. bias:", round(est_mc_plain$rel.abs.bias,3)))
  
  # bootstrap variance estimation
  df_mc_plain <- widetoLong(Y= noisy_mat, mask = mask, X = NULL)
  est_mc_plain$boot_var <- clustered_bootstrap(current_data_realized_long=df_mc_plain, estimator="mc_plain", N=nrow(noisy_mat), T=ncol(noisy_mat), B = 999, est_weights = FALSE)
  print(paste("MC-NNM (Plain) variance:", round(est_mc_plain$boot_var,3)))
  
  est_mc_plain$cp <- CI_test(est_coefficent=est_mc_plain$att.bar, real_coefficent=att.true, est_var=est_mc_plain$boot_var)
  print(paste("MC-NNM (Plain) CP:", round(est_mc_plain$cp,3)))
  
  ## ------ ------ ------ ------ ------
  ## MC-NNM weights (no covariate)
  ## ------ ------ ------ ------ ------
  
  est_mc_weights <- mcnnm(M = obs_mat, mask = mask, W = weights, to_estimate_u = 1, to_estimate_v = 1, lambda_L = 0.0014828, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
  est_mc_weights$Mhat <- est_mc_weights$L + replicate(T,est_mc_weights$u) + t(replicate(N,est_mc_weights$v))
  est_mc_weights$tau <- (est_mc_weights$Mhat-noisy_mat) # estimated treatment effect, Y(AT) - Y(ST)
  est_mc_weights$err <- (est_mc_weights$tau  - shift) # error (wrt to ground truth)
  
  est_mc_weights$msk_err <- est_mc_weights$err*(1-mask) # masked error (wrt to ground truth)
  est_mc_weights$RMSE <- sqrt((1/sum(1-mask)) * sum(est_mc_weights$msk_err^2)) # RMSE on test set (wrt to ground truth)
  print(paste("MC-NNM (weights) RMSE:", round(est_mc_weights$RMSE,3)))
  
  est_mc_weights$att <- apply(est_mc_weights$tau*(1-mask),1,nzmean)[ST]
  est_mc_weights$att.bar <- mean(est_mc_weights$att)
  est_mc_weights$abs.bias <- abs(est_mc_weights$att.bar-att.true)
  est_mc_weights$rel.abs.bias <- est_mc_weights$abs.bias/att.true
  print(paste("MC-NNM (weights) rel. abs. bias:", round(est_mc_weights$rel.abs.bias,3)))
  
  # bootstrap variance estimation
  df_mc_weights <- widetoLong(Y= noisy_mat, mask = mask, X = X)
  est_mc_weights$boot_var <- clustered_bootstrap(current_data_realized_long=df_mc_weights, estimator="mc_weights", N=nrow(noisy_mat), T=ncol(noisy_mat), B = 999, est_weights = TRUE)
  print(paste("MC-NNM (weights) variance:", round(est_mc_weights$boot_var,3)))
  
  est_mc_weights$cp <- CI_test(est_coefficent=est_mc_weights$att.bar, real_coefficent=att.true, est_var=est_mc_weights$boot_var)
  print(paste("MC-NNM (weights) CP:", round(est_mc_weights$cp,3)))
  
  ## ------ ------ ------ ------ ------
  ## MC-NNM + weights + covariate
  ## ------ ------ ------ ------ ------
  
  est_mc_weights_covars <- mcnnm_wc(M = obs_mat, C = X, mask = mask, W = weights, to_estimate_u = 1, to_estimate_v = 1, lambda_L = 0.00108901, lambda_B = 0.100334, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
  est_mc_weights_covars$Mhat <- est_mc_weights_covars$L + X.hat%*%replicate(T,as.vector(est_mc_weights_covars$B)) + replicate(T,est_mc_weights_covars$u) + t(replicate(N,est_mc_weights_covars$v))
  est_mc_weights_covars$tau <- (est_mc_weights_covars$Mhat-noisy_mat) # estimated treatment effect
  est_mc_weights_covars$err <- (est_mc_weights_covars$tau  - shift) # error (wrt to ground truth)
  
  est_mc_weights_covars$msk_err <- est_mc_weights_covars$err*(1-mask) # masked error (wrt to ground truth)
  est_mc_weights_covars$RMSE <- sqrt((1/sum(1-mask)) * sum(est_mc_weights_covars$msk_err^2)) # RMSE on test set (wrt to ground truth)
  print(paste("MC-NNM (weights+covars) RMSE:", round(est_mc_weights_covars$RMSE,3)))
  
  est_mc_weights_covars$att <- apply(est_mc_weights_covars$tau*(1-mask),1,nzmean)[ST]
  est_mc_weights_covars$att.bar <- mean(est_mc_weights_covars$att)
  est_mc_weights_covars$abs.bias <- abs(est_mc_weights_covars$att.bar-att.true)
  est_mc_weights_covars$rel.abs.bias <- est_mc_weights_covars$abs.bias/att.true
  print(paste("MC-NNM (weights+covars) rel. abs. bias:", round(est_mc_weights_covars$rel.abs.bias,3)))
  
  # bootstrap variance estimation
  df_mc_weights_covars <- widetoLong(Y= noisy_mat, mask = mask, X = X)
  est_mc_weights_covars$boot_var <- clustered_bootstrap(current_data_realized_long=df_mc_weights_covars, estimator="mc_weights_covars",N=nrow(noisy_mat), T=ncol(noisy_mat), B = 999, est_weights = TRUE)
  print(paste("MC-NNM (weights + covars) variance:", round(est_mc_weights_covars$boot_var,3)))
  
  est_mc_weights_covars$cp <- CI_test(est_coefficent=est_mc_weights_covars$att.bar, real_coefficent=att.true, est_var=est_mc_weights_covars$boot_var)
  print(paste("MC-NNM (weights + covars) CP", round(est_mc_weights_covars$cp,3)))
  
  ## -----
  ## ADH
  ## -----
  est_model_ADH <- list()
  est_model_ADH$Mhat <- adh_mp_rows(obs_mat, mask)
  est_model_ADH$tau <- (est_model_ADH$Mhat-noisy_mat) # estimated treatment effect
  est_model_ADH$err <- (est_model_ADH$tau  - shift) # error (wrt to ground truth)
  
  est_model_ADH$msk_err <- est_model_ADH$err*(1-mask) # masked error (wrt to ground truth)
  est_model_ADH$RMSE <- sqrt((1/sum(1-mask)) * sum(est_model_ADH$msk_err^2)) # RMSE on test set (wrt to ground truth)
  print(paste("ADH RMSE:", round(est_model_ADH$RMSE,3)))
  
  est_model_ADH$att <- apply(est_model_ADH$tau*(1-mask),1,nzmean)[ST]
  est_model_ADH$att.bar <- mean(est_model_ADH$att)
  est_model_ADH$abs.bias <- abs(est_model_ADH$att.bar-att.true)
  est_model_ADH$rel.abs.bias <- est_model_ADH$abs.bias/att.true
  print(paste("ADH rel. abs. bias:", round(est_model_ADH$rel.abs.bias,3)))
  
  # bootstrap variance estimation
  df_ADH <- widetoLong(Y= noisy_mat, mask = mask, X = NULL)
  est_model_ADH$boot_var <- clustered_bootstrap(current_data_realized_long=df_ADH, estimator="ADH", N=nrow(noisy_mat), T=ncol(noisy_mat), B = 999, est_weights = FALSE)
  print(paste("ADH variance:", round(est_model_ADH$boot_var,3)))
  
  est_model_ADH$cp <- CI_test(est_coefficent=est_model_ADH$att.bar, real_coefficent=att.true, est_var=est_model_ADH$boot_var)
  print(paste("ADH CP:", round(est_model_ADH$cp,3)))
  
  ## -----
  ## DID
  ## -----
  est_model_DID <- list()
  est_model_DID$Mhat <- DID(obs_mat, mask)
  est_model_DID$tau <- (est_model_DID$Mhat-noisy_mat) # estimated treatment effect
  est_model_DID$err <- (est_model_DID$tau  - shift) # error (wrt to ground truth)
  
  est_model_DID$msk_err <- est_model_DID$err*(1-mask) # masked error (wrt to ground truth)
  est_model_DID$RMSE <- sqrt((1/sum(1-mask)) * sum(est_model_DID$msk_err^2)) # RMSE on test set (wrt to ground truth)
  print(paste("DID RMSE:", round(est_model_DID$RMSE,3)))
  
  est_model_DID$att <- apply(est_model_DID$tau*(1-mask),1,nzmean)[ST]
  est_model_DID$att.bar <- mean(est_model_DID$att)
  est_model_DID$abs.bias <- abs(est_model_DID$att.bar-att.true)
  est_model_DID$rel.abs.bias <- est_model_DID$abs.bias/att.true
  print(paste("DID rel. abs. bias:", round(est_model_DID$rel.abs.bias,3)))
  
  # bootstrap variance estimation
  df_DID <- widetoLong(Y= noisy_mat, mask = mask, X = NULL)
  est_model_DID$boot_var <- clustered_bootstrap(current_data_realized_long=df_DID, estimator="DID", N=nrow(noisy_mat), T=ncol(noisy_mat), B = 999, est_weights = FALSE)
  print(paste("DID variance:", round(est_model_DID$boot_var,3)))
  
  est_model_DID$cp <- CI_test(est_coefficent=est_model_DID$att.bar, real_coefficent=att.true, est_var=est_model_DID$boot_var)
  print(paste("DID CP:", round(est_model_DID$cp,3)))
  
  ## ---------------
  ## IFEs
  ## ---------------
  
  est_model_IFE <- list()
  est_model_IFE$Mhat <- IFE(obs_mat, mask, k=2)
  est_model_IFE$tau <- (est_model_IFE$Mhat-noisy_mat) # estimated treatment effect
  est_model_IFE$err <- (est_model_IFE$tau  - shift) # error (wrt to ground truth)
  
  est_model_IFE$msk_err <- est_model_IFE$err*(1-mask) # masked error (wrt to ground truth)
  est_model_IFE$RMSE <- sqrt((1/sum(1-mask)) * sum(est_model_IFE$msk_err^2)) # RMSE on test set (wrt to ground truth)
  print(paste("IFE RMSE:", round(est_model_IFE$RMSE,3)))
  
  est_model_IFE$att <- apply(est_model_IFE$tau*(1-mask),1,nzmean)[ST]
  est_model_IFE$att.bar <- mean(est_model_IFE$att)
  est_model_IFE$abs.bias <- abs(est_model_IFE$att.bar-att.true)
  est_model_IFE$rel.abs.bias <- est_model_IFE$abs.bias/att.true
  print(paste("IFE rel. abs. bias:", round(est_model_IFE$rel.abs.bias,3)))
  
  # bootstrap variance estimation
  df_IFE <- widetoLong(Y= noisy_mat, mask = mask, X = NULL)
  est_model_IFE$boot_var <- clustered_bootstrap(current_data_realized_long=df_IFE, estimator="IFE", N=nrow(noisy_mat), T=ncol(noisy_mat), B = 999, est_weights = FALSE)
  print(paste("IFE variance:", round(est_model_IFE$boot_var,3)))
  
  est_model_IFE$cp <- CI_test(est_coefficent=est_model_IFE$att.bar, real_coefficent=att.true, est_var=est_model_IFE$boot_var)
  print(paste("IFE CP:", round(est_model_IFE$cp,3)))
  
  # cleanup
  rm(A,B,C,D,e,mask,noise,noisy_mat,obs_mat,true_mat,shifted_mat,sigma_mat_N,sigma_mat_N_treat,sigma_mat_R,sigma_mat_R_treat,treat_mat,X,X.hat,weights,z_weights,df_DID,df_ADH,df_mc_weights_covars,df_mc_plain, df_IFE)
  gc()
  cat(paste("Done with simulation run number",n, "\n"))
  return(list("N"=N, "T"=T, "R"=R, "T0"=T0, "noise_sc"=noise_sc,"delta_sc"=delta_sc, "gamma_sc"=gamma_sc,"beta_sc"=beta_sc, "shift_sc"=shift_sc,"fr_obs"= fr_obs, 
              "est_mc_plain_RMSE"=est_mc_plain$RMSE,"est_mc_plain_abs_bias"=est_mc_plain$abs.bias,"est_mc_plain_rel_abs_bias"=est_mc_plain$rel.abs.bias,"est_mc_plain_cp"=est_mc_plain$cp,"est_mc_plain_boot_var"=est_mc_plain$boot_var,
              "est_mc_weights_RMSE"=est_mc_weights$RMSE,"est_mc_weights_abs_bias"=est_mc_weights$abs.bias,"est_mc_weights_rel_abs_bias"=est_mc_weights$rel.abs.bias,"est_mc_weights_cp"=est_mc_weights$cp,"est_mc_weights_boot_var"=est_mc_weights$boot_var,
              "est_mc_weights_covars_RMSE"=est_mc_weights_covars$RMSE,"est_mc_weights_covars_abs_bias"=est_mc_weights_covars$abs.bias,"est_mc_weights_covars_rel_abs_bias"=est_mc_weights_covars$rel.abs.bias,"est_mc_weights_covars_cp"=est_mc_weights_covars$cp,"est_mc_weights_covars_boot_var"=est_mc_weights_covars$boot_var,
              "est_model_ADH_RMSE"=est_model_ADH$RMSE,"est_model_ADH_abs_bias"=est_model_ADH$abs.bias,"est_model_ADH_rel_abs_bias"=est_model_ADH$rel.abs.bias,"est_model_ADH_cp"=est_model_ADH$cp,"est_model_ADH_boot_var"=est_model_ADH$boot_var,
              "est_model_DID_RMSE"=est_model_DID$RMSE,"est_model_DID_abs_bias"=est_model_DID$abs.bias,"est_model_DID_rel_abs_bias"=est_model_DID$rel.abs.bias,"est_model_DID_cp"=est_model_DID$cp,"est_model_DID_boot_var"=est_model_DID$boot_var,
              "est_model_IFE_RMSE"=est_model_IFE$RMSE,"est_model_IFE_abs_bias"=est_model_IFE$abs.bias,"est_model_IFE_rel_abs_bias"=est_model_IFE$rel.abs.bias,"est_model_IFE_cp"=est_model_IFE$cp,"est_model_IFE_boot_var"=est_model_IFE$boot_var))
}

# define settings for simulation
settings <- expand.grid("NT"=c(40**2),
                        "N_t" =c(0.4,0.6,0.8),
                        "R"=c(10,20,30,40))

args <- commandArgs(trailingOnly = TRUE) # command line arguments
thisrun <- settings[as.numeric(args[1]),] 

N <- sqrt(as.numeric(thisrun[1])) # Number of units
T <- sqrt(as.numeric(thisrun[1]))  # Number of time-periods
T0 <- ceiling(T*0.5)  # initial time period in staggered adoption setting
N_t <- ceiling(N*as.numeric(thisrun[2]))  # Number of ST units
R <- as.numeric(thisrun[3])

noise_sc <- 0.01 # Noise scale 
delta_sc <- 0.1 # delta scale
gamma_sc <- 0.1 # gamma scale
beta_sc <- 0.1 # beta scale
shift_sc <- 0.01 # shift scale

n.runs <- 1000 # Num. simulation runs

output_dir <- './outputs/'
simulation_version <- paste0(format(Sys.time(), "%Y%m%d"),"/")
if(!dir.exists(output_dir)){
  print(paste0('create folder for outputs at: ', output_dir))
  dir.create(output_dir)
}
output_dir <- paste0(output_dir, simulation_version)
if(!dir.exists(output_dir)){
  print(paste0('create folder for outputs at: ', output_dir))
  dir.create(output_dir)
}

setting <- paste0("N = ", N, ", T = ", T, ", R = ",R, "T0 = " ,T0, "N_t", N_t,", noise_sc = ",noise_sc, ", delta_sc = ",delta_sc, ", gamma_sc = ",gamma_sc, ", beta_sc = ", beta_sc, ", shift_sc = ", shift_sc)
tic(print(paste0("setting: ",setting)))

results <- foreach(i = 1:n.runs, .combine='rbind', .packages =c("MCPanel","matrixStats","Matrix","MASS","data.table","reshape","reshape2","emfactor"), .verbose = FALSE) %dopar% {
  MCsim(N,T,R,T0,N_t,noise_sc,delta_sc,gamma_sc,beta_sc,shift_sc,n=i)
}
results
saveRDS(results, paste0(output_dir,"results_","N_",N,"_T_",T,"_R_", R,"_T0_",T0, "_N_t", N_t, "_noise_sc_",noise_sc,"_shift_sc_",shift_sc,"_n_",n.runs,".rds"))

print(toc())

if(doMPI){
  closeCluster(cl) # close down MPIcluster
  mpi.finalize()
}else{
  stopCluster(cl)
}