##################################
# Matrix completion simulations  #
##################################

library(MCPanel)
library(matrixStats)
library(boot)
library(Matrix)

library(ggplot2)
library(latex2exp)
library(dplyr)

# Setup parallel processing
library(parallel)
library(doParallel)

cores <- detectCores()

cl <- parallel::makeForkCluster(cores)

doParallel::registerDoParallel(cores) # register cores (<p)

RNGkind("L'Ecuyer-CMRG") # ensure random number generation

## Set random seed
set.seed(10)

MCsim <- function(N,T,R,noise_sc,delta_sc,gamma_sc,beta_sc,effect_size){
  
  setting <- paste0("N = ", N, ", T = ", T, ", R = ",R, ", noise_sc = ",noise_sc, ", delta_sc = ",delta_sc, ", gamma_sc = ",gamma_sc, ", beta_sc = ", beta_sc, ", effect_size = ", effect_size)
  print(setting)
  
  # Create Matrices
  A <- replicate(R,rnorm(N))
  B <- replicate(T,rnorm(R))
  X <- replicate(T,rnorm(N))
  delta <- delta_sc*rnorm(N)
  gamma <- gamma_sc*rnorm(T)
  beta <- beta_sc*rnorm(T)
  noise <- noise_sc*replicate(T,rnorm(N))
  
  # True outcome model
  true_mat_0 <- A %*% B + X%*%replicate(T,as.vector(beta)) + replicate(T,delta) + t(replicate(N,gamma)) # potential outcome under control
  true_mat_1 <- true_mat_0 + effect_size # potential outcome under treatment
  
  noisy_mat <- true_mat_1 + noise # we want to estimate p.o. under treatment in retrospective analysis

  # Assignment mechanism 
  C <- replicate((N-1),rbinom(N,1,delta_sc))
  D <- replicate(T,rbinom((T-1),1,gamma_sc)) 
  
  e <-plogis(C%*%D + X%*%replicate(T,as.vector(beta)) + replicate(T,delta) + t(replicate(N,gamma)) + noise) #incl. noise
  treat_mat <- matrix(rbinom(N*T,1,e),N,T)
  
  # block structure (fill forwards)
  for (i in 1:N){ 
    for (j in 1:(T-1)) {
      if (treat_mat[i,j]==1) {treat_mat[i,(j+1)]=1}
      else {treat_mat[i,j]=treat_mat[i,j]}
    }
  }

  # Observed outcome

  mask <- treat_mat # 0s missing and to be imputed; 1s are observed
  
  obs_mat <- noisy_mat * mask 
  
  # Impute endogenous covariate values
  
  est_model_endogenous <- mcnnm_cv(M = X, mask = mask, W = matrix(1, nrow(mask),ncol(mask)), to_estimate_u = 1, to_estimate_v = 1, 
                                   num_lam_L = 30, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 5, is_quiet = 1) # outcome is covariate, prop. weights are equal
  est_model_endogenous$Mhat <- est_model_endogenous$L + replicate(T,est_model_endogenous$u) + t(replicate(N,est_model_endogenous$v))

  X.hat <- X*mask + est_model_endogenous$Mhat*(1-mask) # only endogenous values imputed
  
  # Estimate propensity weights by matrix completion

  est_model_pweights <- mcnnm_wc_cv(M = (1-mask), C = X.hat, mask =  matrix(1, nrow(mask),ncol(mask)), # no missing entries
                                    W = matrix(1, nrow(mask),ncol(mask)), to_normalize = 1, to_estimate_u = 1, to_estimate_v = 1, num_lam_L = 30, num_lam_B = 30, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 5, is_quiet = 1) # use X with imputed endogenous values
  est_model_pweights$Mhat <- est_model_pweights$L + X.hat%*%replicate(T,as.vector(est_model_pweights$B)) + replicate(T,est_model_pweights$u) + t(replicate(N,est_model_pweights$v)) # use X with imputed endogenous values
  
  source('boundProbs.R')
  
  weights <- (1-boundProbs(est_model_pweights$Mhat))/boundProbs(est_model_pweights$Mhat)
  z_weights <- c(rev(SSlogis(1:t0, Asym = 1, xmid = 0.85, scal = 8)),
                 SSlogis(1:(ncol(mask)-t0), Asym = 1, xmid = 0.85, scal = 8))
  weights <- weights%*%diag(z_weights) # elapsed-time weighting

  # Function for bootstrap
  BootTraj <- function(impact,indices,t0,treat_indices,start=1) {
    att <- rowMeans(impact[,start:(t0-1)])[indices] 
    return(mean(att[treat_indices])) # avg over treated units
  }
  
  ## -----
  ## ADH
  ## -----
  est_model_ADH <- list()
  est_model_ADH$Mhat <- adh_mp_rows(obs_mat, mask, niter=200, rel_tol = 0.001)
  est_model_ADH$impact <- (est_model_ADH$Mhat-noisy_mat) # estimated treatment effect
  est_model_ADH$err <- (est_model_ADH$Mhat - true_mat_1) # error (wrt to ground truth)

  est_model_ADH$msk_err <- est_model_ADH$err*(1-mask) # masked error (wrt to ground truth)
  est_model_ADH$test_RMSE <- sqrt((1/sum(1-mask)) * sum(est_model_ADH$msk_err^2)) # RMSE on test set (wrt to ground truth)
  print(paste("DID RMSE:", round(est_model_ADH$test_RMSE,3)))
  
  est_model_ADH$boot <- boot(est_model_ADH$impact, 
                            BootTraj, 
                            t0=t0,
                            treat_indices=treat_indices,
                            R=999,
                            parallel = "multicore") 
  
  est_model_ADH$boot_t0 <- est_model_ADH$boot$t0
  est_model_ADH$boot_ci <- boot.ci(est_model_ADH$boot ,type=c("perc"))$percent[-c(1:3)] 
  est_model_ADH$bias <- est_model_ADH$boot_t0-effect_size
  est_model_ADH$cp <- as.numeric((est_model_ADH$boot_ci[1] < effect_size) & (est_model_ADH$boot_ci[2] > effect_size))
  
  ## -----
  ## DID
  ## -----
  est_model_DID <- list()
  est_model_DID$Mhat <- t(DID(t(obs_mat), t(mask)))
  est_model_DID$impact <- (est_model_DID$Mhat-noisy_mat) # estimated treatment effect
  est_model_DID$err <- (est_model_DID$Mhat - true_mat_1) # error (wrt to ground truth)

  est_model_DID$msk_err <- est_model_DID$err*(1-mask) # masked error (wrt to ground truth)
  est_model_DID$test_RMSE <- sqrt((1/sum(1-mask)) * sum(est_model_DID$msk_err^2)) # RMSE on test set (wrt to ground truth)
  print(paste("DID RMSE:", round(est_model_DID$test_RMSE,3)))
  
  est_model_DID$boot <- boot(est_model_DID$impact, 
                             BootTraj, 
                             t0=t0,
                             treat_indices=treat_indices,
                             R=999,
                             parallel = "multicore") 
  
  est_model_DID$boot_t0 <- est_model_DID$boot$t0
  est_model_DID$boot_ci <- boot.ci(est_model_DID$boot ,type=c("perc"))$percent[-c(1:3)] 
  est_model_DID$bias <- est_model_DID$boot_t0-effect_size
  est_model_DID$cp <- as.numeric((est_model_DID$boot_ci[1] < effect_size) & (est_model_DID$boot_ci[2] > effect_size))
  
  ## ------ ------ ------ ------ ------
  ## MC-NNM plain (no weighting, no covariates)
  ## ------ ------ ------ ------ ------
  
  est_mc_plain <- mcnnm_cv(M = obs_mat, mask = mask, W = matrix(1, nrow(mask),ncol(mask)), to_estimate_u = 1, to_estimate_v = 1, num_lam_L = 30, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 5, is_quiet = 1)
  est_mc_plain$Mhat <- est_mc_plain$L + replicate(T,est_mc_plain$u) + t(replicate(N,est_mc_plain$v))
  est_mc_plain$impact <- (est_mc_plain$Mhat-noisy_mat) # estimated treatment effect
  est_mc_plain$err <- (est_mc_plain$Mhat - true_mat_1) # error (wrt to ground truth)

  est_mc_plain$msk_err <- est_mc_plain$err*(1-mask) # masked error (wrt to ground truth)
  est_mc_plain$test_RMSE <- sqrt((1/sum(1-mask)) * sum(est_mc_plain$msk_err^2)) # RMSE on test set (wrt to ground truth)
  print(paste("MC-NNM (Plain) RMSE:", round(est_mc_plain$test_RMSE,3)))

  est_mc_plain$boot <- boot(est_mc_plain$impact, 
                        BootTraj, 
                        t0=t0,
                        treat_indices=treat_indices,
                        R=999,
                        parallel = "multicore") 
  
  est_mc_plain$boot_t0 <- est_mc_plain$boot$t0
  est_mc_plain$boot_ci <- boot.ci(est_mc_plain$boot ,type=c("perc"))$percent[-c(1:3)] 
  est_mc_plain$bias <- est_mc_plain$boot_t0-effect_size
  est_mc_plain$cp <- as.numeric((est_mc_plain$boot_ci[1] < effect_size) & (est_mc_plain$boot_ci[2] > effect_size))
  
  ## ------ ------ ------ ------ ------
  ## MC-NNM + weights (no covariates)
  ## ------ ------ ------ ------ ------
  
  est_mc_weights <- mcnnm_cv(M = obs_mat, mask = mask, W = weights, to_estimate_u = 1, to_estimate_v = 1, num_lam_L = 30, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 5, is_quiet = 1)
  est_mc_weights$Mhat <- est_mc_weights$L + replicate(T,est_mc_weights$u) + t(replicate(N,est_mc_weights$v))
  est_mc_weights$impact <- (est_mc_weights$Mhat-noisy_mat) # estimated treatment effect
  est_mc_weights$err <- (est_mc_weights$Mhat - true_mat_1) # error (wrt to ground truth)

  est_mc_weights$msk_err <- est_mc_weights$err*(1-mask) # masked error (wrt to ground truth)
  est_mc_weights$test_RMSE <- sqrt((1/sum(1-mask)) * sum(est_mc_weights$msk_err^2)) # RMSE on test set (wrt to ground truth)
  print(paste("MC-NNM (weights) RMSE:", round(est_mc_weights$test_RMSE,3)))
  
  est_mc_weights$boot <- boot(est_mc_weights$impact, 
                            BootTraj, 
                            t0=t0,
                            treat_indices=treat_indices,
                            R=999,
                            parallel = "multicore") 
  
  est_mc_weights$boot_t0 <- est_mc_weights$boot$t0
  est_mc_weights$boot_ci <- boot.ci(est_mc_weights$boot ,type=c("perc"))$percent[-c(1:3)] 
  est_mc_weights$bias <- est_mc_weights$boot_t0-effect_size
  est_mc_weights$cp <- as.numeric((est_mc_weights$boot_ci[1] < effect_size) & (est_mc_weights$boot_ci[2] > effect_size))
  
  ## ------ ------ ------ ------ ------
  ## MC-NNM + weights + covariates
  ## ------ ------ ------ ------ ------
  
  est_mc_weights_covars <- mcnnm_wc_cv(M = obs_mat, C = X, mask = mask, W = weights, to_estimate_u = 1, to_estimate_v = 1, num_lam_L = 30, num_lam_B = 30, niter = 1000, rel_tol = 1e-05, cv_ratio = 0.8, num_folds = 5, is_quiet = 1)
  est_mc_weights_covars$Mhat <- est_mc_weights_covars$L + X.hat%*%replicate(T,as.vector(est_mc_weights_covars$B)) + replicate(T,est_mc_weights_covars$u) + t(replicate(N,est_mc_weights_covars$v))
  est_mc_weights_covars$impact <- (est_mc_weights_covars$Mhat-noisy_mat) # estimated treatment effect
  est_mc_weights_covars$err <- (est_mc_weights_covars$Mhat - true_mat_1) # error (wrt to ground truth)

  est_mc_weights_covars$msk_err <- est_mc_weights_covars$err*(1-mask) # masked error (wrt to ground truth)
  est_mc_weights_covars$test_RMSE <- sqrt((1/sum(1-mask)) * sum(est_mc_weights_covars$msk_err^2)) # RMSE on test set (wrt to ground truth)
  print(paste("MC-NNM (covariates) RMSE:", round(est_mc_weights_covars$test_RMSE,3)))
  
  est_mc_weights_covars$boot <- boot(est_mc_weights_covars$impact, 
                              BootTraj, 
                              t0=t0,
                              treat_indices=treat_indices,
                              R=999,
                              parallel = "multicore") 
  
  est_mc_weights_covars$boot_t0 <- est_mc_weights_covars$boot$t0
  est_mc_weights_covars$boot_ci <- boot.ci(est_mc_weights_covars$boot ,type=c("perc"))$percent[-c(1:3)] 
  est_mc_weights_covars$bias <- est_mc_weights_covars$boot_t0-effect_size
  est_mc_weights_covars$cp <- as.numeric((est_mc_weights_covars$boot_ci[1] < effect_size) & (est_mc_weights_covars$boot_ci[2] > effect_size))
    
  return(list("setting"=setting, 
              "est_mc_plain_RMSE"=est_mc_plain$test_RMSE,"est_mc_plain_bias"=est_mc_plain$bias,"est_mc_plain_cp"=est_mc_plain$cp,
              "est_mc_weights_RMSE"=est_mc_weights$test_RMSE,"est_mc_weights_bias"=est_mc_weights$bias,"est_mc_weights_cp"=est_mc_weights$cp,
              "est_mc_weights_covars_RMSE"=est_mc_weights_covars$test_RMSE,"est_mc_weights_covars_bias"=est_mc_weights_covars$bias,"est_mc_weights_covars_cp"=est_mc_weights_covars$cp,
              "est_model_ADH_RMSE"=est_model_ADH$test_RMSE,"est_model_ADH_bias"=est_model_ADH$bias,"est_model_ADH_cp"=est_model_ADH$cp,
              "est_model_DID_RMSE"=est_model_DID$test_RMSE,"est_model_DID_bias"=est_model_DID$bias,"est_model_DID_cp"=est_model_DID$cp))
}

# set the seed
set.seed(10)

# define settings for simulation
settings <- expand.grid("N"=c(50,100,500),
                        "T"=c(50,100,500),
                        "noise_sc"=c(0.001,0.01,0.1),
                        "effect_size"=c(0.001,0.01,0.1))

settings$R <- NA
settings[settings$N==50 | settings$T==50,]$R <- 15 
settings[settings$N==100 | settings$T==100,]$R <- 25 
settings[settings$N==500 & settings$T==500,]$R <- 50
  
args <- as.numeric(commandArgs(trailingOnly = TRUE)) # command line arguments
thisrun <- settings[args,] 

N <- as.numeric(thisrun[1]) # Number of units
T <- as.numeric(thisrun[2]) # Number of time-periods

noise_sc <- as.numeric(thisrun[3]) # Noise scale 
delta_sc <- 0.01 # delta scale
gamma_sc <- 0.05 # gamma scale
beta_sc <- 0.1 # beta scale
effect_size <- as.numeric(thisrun[4])
R <- as.numeric(thisrun[5])

number_T0 <- 4
T0 <- ceiling(T*((1:number_T0)*2-1)/(2*number_T0))

n <- 1000 # Num. simulation runs

results <- foreach(t0 =T0, .combine='rbind') %dopar% {
  replicate(n,MCsim(N,T,R,noise_sc,delta_sc,gamma_sc,beta_sc,effect_size))
}
results <- matrix(unlist(results), ncol = n, byrow = FALSE) # coerce into matrix
saveRDS(results, paste0("results_","N_",N,"_T_",T,"_R_", R,"_noise_sc_",noise_sc,"_effect_size_",effect_size,"_n_",n,".rds"))