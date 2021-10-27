##################################
# Matrix completion simulations  #
##################################

library(MCPanel)
library(matrixStats)
library(boot)
library(Matrix)
library(tictoc)

source('boundProbs.R')
source('bootTraj.R')

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

MCsim <- function(N,T,R,noise_sc,delta_sc,gamma_sc,beta_sc,effect_size,n){
  
  # set the seed
  print(paste0("run number: ", n))
  set.seed(n, "L'Ecuyer-CMRG")
  
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
  
  C <- replicate(R,rbinom(N,1,0.3))
  D <- replicate(T,rbinom(R,1,0.3))
  
  e <-1/(1+exp(C%*%D + X%*%replicate(T,as.vector(beta)) + replicate(T,delta) + t(replicate(N,gamma)) + noise)) #incl. noise
  
  treat_mat <- matrix(rbinom(N*T,1,e),N,T) # 0s missing and to be imputed (LT); 1s are observed (AT) 
  while(any(rowSums(treat_mat)<2) || max(rowSums(treat_mat))<T){ # generate new treat_mat to ensure all units treated for at least 2 periods and that there are AT units
    treat_mat <- matrix(rbinom(N*T,1,e),N,T)  
    
    # block structure (fill forwards)
    for (i in 1:N){ 
      for (j in 1:(T-1)) {
        if (treat_mat[i,j]==1) {treat_mat[i,(j+1)]=1}
        else {treat_mat[i,j]=treat_mat[i,j]}
      }
    }
  }
  
  fr_obs <- sum(treat_mat)/(N*T) # store fraction observed entries
  
  # Get vector of initial treatment periods
  
  T0 <- aggregate(col ~ row,
                  data = which(treat_mat == 1, arr.ind = T),
                  FUN = function(x) x[1])$col
  
  # Observed outcome

  mask <- treat_mat 
  obs_mat <- noisy_mat * mask 
  
  # CV lambda values
  
  lambda_L_endog <- as.numeric(sort(c(0.0201769,0.0202961,0.0258464,0.0239927,0.0202191,0.0216138,0.021132,0.0216067,0.0232151,0.021332,0.0201098,0.0203622,0.0248053,0.0219271,0.0196738,0.0237891)))
  lambda_L_pweights <- as.numeric(sort(c(0,8.69584e-05,9.64993e-05,0.000886265,8.71857e-05,9.21348e-05,0.000839902,9.90469e-05,0.000888573,0.000879063,9.11552e-05,8.59036e-05,9.29274e-05,9.0942e-05,9.83006e-05)))
  lambda_B_pweights <- as.numeric(sort(c(0,0.00886265,0.00964993,0.000871857,0.00921348,8.39902e-06,0.000990469,0.00888573,0.000860964,0.000859036,0.00929274,0.0090942,0.00983006)))
  lambda_L_plain <- as.numeric(sort(c(0,0.00149916,0.00153424,0.0090925,0.00344375,0.000557275,0.000584998,0.00331317,0.00146644,0.00391252,0.00407219,0.00133857,0.0038706,0.00133738,0.000108937,0.00141742)))
  lambda_L_weights <- as.numeric(sort(c(0.00153424,0.00355506,0.0090925,0.00344375,0.00132151,0.00328972,0.00331313,0.00347747,0.00391252,0.0040722,0.00317424,0.00142967,0.00917863,0.00317143,0.0014527,0.00336123)))
  lambda_L_covars <- as.numeric(sort(c(0,0.00115052,0.00112421,0.012125,0.00108901,0.0010477,0.00109967,0.000123725,0.00128774,0.00100378, 0.00122399,0.0010029,0.00106291)))
  lambda_B_covars <- as.numeric(sort(c(0.0115052,0.112421,0.012125,0.108901,0.00990989,0.010403,0.10477,0.109967,0.0123725,0.128774,0.100378,0.010721,0.122399,0.10029,0.0108937,0.106291)))
  
  # Impute endogenous covariate values
  
  est_model_endogenous <- mcnnm(M = X, mask = mask, W = matrix(1, nrow(mask),ncol(mask)), to_estimate_u = 1, to_estimate_v = 1, 
                                lambda_L = median(lambda_L_endog), niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]] # outcome is covariate, prop. weights are equal
  est_model_endogenous$Mhat <- est_model_endogenous$L + replicate(T,est_model_endogenous$u) + t(replicate(N,est_model_endogenous$v))

  X.hat <- X*mask + est_model_endogenous$Mhat*(1-mask) # only endogenous values imputed
  
  # Estimate propensity weights by matrix completion

  est_model_pweights <- mcnnm_wc(M = (1-mask), C = X.hat, mask =  matrix(1, nrow(mask),ncol(mask)), # no missing entries
                                    W = matrix(1, nrow(mask),ncol(mask)), to_normalize = 1, to_estimate_u = 1, to_estimate_v = 1, lambda_L = median(lambda_L_pweights), lambda_B = median(lambda_B_pweights), niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]] # use X with imputed endogenous values
  est_model_pweights$Mhat <- est_model_pweights$L + X.hat%*%replicate(T,as.vector(est_model_pweights$B)) + replicate(T,est_model_pweights$u) + t(replicate(N,est_model_pweights$v)) # use X with imputed endogenous values
  
  weights <- (1-boundProbs(est_model_pweights$Mhat))/boundProbs(est_model_pweights$Mhat)
  
  z_weights <- matrix(0,N,T,byrow = TRUE)
  
  for(i in 1:N){
    z_weights[i,] <- c(plogis(T0[i]:1, scale=8),plogis(1:(ncol(mask)-T0[i]),scale=8))
  }
         
  weights <- (1-diag(z_weights)*weights)/(diag(z_weights)*weights) # elapsed-time weighting
  
  ## ------ ------ ------ ------ ------
  ## MC-NNM plain (no weighting, no covariates)
  ## ------ ------ ------ ------ ------
  
  est_mc_plain <- mcnnm(M = obs_mat, mask = mask, W = matrix(1, nrow(mask),ncol(mask)), to_estimate_u = 1, to_estimate_v = 1, lambda_L = median(lambda_L_plain), niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
  est_mc_plain$Mhat <- est_mc_plain$L + replicate(T,est_mc_plain$u) + t(replicate(N,est_mc_plain$v))
  est_mc_plain$impact <- (est_mc_plain$Mhat-noisy_mat) # estimated treatment effect
  est_mc_plain$err <- (est_mc_plain$Mhat - true_mat_1) # error (wrt to ground truth)

  est_mc_plain$msk_err <- est_mc_plain$err*(1-mask) # masked error (wrt to ground truth)
  est_mc_plain$test_RMSE <- sqrt((1/sum(1-mask)) * sum(est_mc_plain$msk_err^2)) # RMSE on test set (wrt to ground truth)
  print(paste("MC-NNM (Plain) RMSE:", round(est_mc_plain$test_RMSE,3)))

  est_mc_plain$boot <- boot(est_mc_plain$impact, 
                        BootTraj, 
                        mask=mask,
                        R=999,
                        parallel = "multicore") 
  
  est_mc_plain$boot_t0 <- est_mc_plain$boot$t0
  est_mc_plain$boot_ci <- boot.ci(est_mc_plain$boot ,type=c("perc"))$percent[-c(1:3)] 
  est_mc_plain$bias <- est_mc_plain$boot_t0-effect_size
  est_mc_plain$cp <- as.numeric((est_mc_plain$boot_ci[1] < effect_size) & (est_mc_plain$boot_ci[2] > effect_size))
  
  ## ------ ------ ------ ------ ------
  ## MC-NNM + weights (no covariates)
  ## ------ ------ ------ ------ ------
  
  est_mc_weights <- mcnnm(M = obs_mat, mask = mask, W = weights, to_estimate_u = 1, to_estimate_v = 1, lambda_L = median(lambda_L_weights), niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
  est_mc_weights$Mhat <- est_mc_weights$L + replicate(T,est_mc_weights$u) + t(replicate(N,est_mc_weights$v))
  est_mc_weights$impact <- (est_mc_weights$Mhat-noisy_mat) # estimated treatment effect
  est_mc_weights$err <- (est_mc_weights$Mhat - true_mat_1) # error (wrt to ground truth)

  est_mc_weights$msk_err <- est_mc_weights$err*(1-mask) # masked error (wrt to ground truth)
  est_mc_weights$test_RMSE <- sqrt((1/sum(1-mask)) * sum(est_mc_weights$msk_err^2)) # RMSE on test set (wrt to ground truth)
  print(paste("MC-NNM (weights) RMSE:", round(est_mc_weights$test_RMSE,3)))
  
  est_mc_weights$boot <- boot(est_mc_weights$impact, 
                            BootTraj, 
                            mask=mask,
                            R=999,
                            parallel = "multicore") 
  
  est_mc_weights$boot_t0 <- est_mc_weights$boot$t0
  est_mc_weights$boot_ci <- boot.ci(est_mc_weights$boot ,type=c("perc"))$percent[-c(1:3)] 
  est_mc_weights$bias <- est_mc_weights$boot_t0-effect_size
  est_mc_weights$cp <- as.numeric((est_mc_weights$boot_ci[1] < effect_size) & (est_mc_weights$boot_ci[2] > effect_size))
  
  ## ------ ------ ------ ------ ------
  ## MC-NNM + weights + covariates
  ## ------ ------ ------ ------ ------
  
  est_mc_weights_covars <- mcnnm_wc(M = obs_mat, C = X, mask = mask, W = weights, to_estimate_u = 1, to_estimate_v = 1, lambda_L = median(lambda_L_covars), lambda_B = median(lambda_B_covars), niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
  est_mc_weights_covars$Mhat <- est_mc_weights_covars$L + X.hat%*%replicate(T,as.vector(est_mc_weights_covars$B)) + replicate(T,est_mc_weights_covars$u) + t(replicate(N,est_mc_weights_covars$v))
  est_mc_weights_covars$impact <- (est_mc_weights_covars$Mhat-noisy_mat) # estimated treatment effect
  est_mc_weights_covars$err <- (est_mc_weights_covars$Mhat - true_mat_1) # error (wrt to ground truth)

  est_mc_weights_covars$msk_err <- est_mc_weights_covars$err*(1-mask) # masked error (wrt to ground truth)
  est_mc_weights_covars$test_RMSE <- sqrt((1/sum(1-mask)) * sum(est_mc_weights_covars$msk_err^2)) # RMSE on test set (wrt to ground truth)
  print(paste("MC-NNM (weights + covariates) RMSE:", round(est_mc_weights_covars$test_RMSE,3)))
  
  est_mc_weights_covars$boot <- boot(est_mc_weights_covars$impact, 
                              BootTraj, 
                              mask=mask,
                              R=999,
                              parallel = "multicore") 
  
  est_mc_weights_covars$boot_t0 <- est_mc_weights_covars$boot$t0
  est_mc_weights_covars$boot_ci <- boot.ci(est_mc_weights_covars$boot ,type=c("perc"))$percent[-c(1:3)] 
  est_mc_weights_covars$bias <- est_mc_weights_covars$boot_t0-effect_size
  est_mc_weights_covars$cp <- as.numeric((est_mc_weights_covars$boot_ci[1] < effect_size) & (est_mc_weights_covars$boot_ci[2] > effect_size))
  
  # cleanup
  rm(A,C,e,mask,noise,noisy_mat,obs_mat,treat_mat,true_mat_0,true_mat_1,X,X.hat,weights,z_weights)
  gc()
  cat(paste("Done with simulation run number",n, "\n"))
    
  return(list("N"=N, "T"=T, "R"=R, "noise_sc"=noise_sc,"delta_sc"=delta_sc, "gamma_sc"=gamma_sc,"beta_sc"=beta_sc, "effect_size"=effect_size, "fr_obs"= fr_obs, 
              "est_mc_plain_RMSE"=est_mc_plain$test_RMSE,"est_mc_plain_bias"=est_mc_plain$bias,"est_mc_plain_cp"=est_mc_plain$cp,
              "est_mc_weights_RMSE"=est_mc_weights$test_RMSE,"est_mc_weights_bias"=est_mc_weights$bias,"est_mc_weights_cp"=est_mc_weights$cp,
              "est_mc_weights_covars_RMSE"=est_mc_weights_covars$test_RMSE,"est_mc_weights_covars_bias"=est_mc_weights_covars$bias,"est_mc_weights_covars_cp"=est_mc_weights_covars$cp))
}

# define settings for simulation
settings <- expand.grid("NT"=c(1600,2500,3600),
                        "noise_sc"=c(0.1,0.2,0.3),
                        "effect_size"=c(0.001,0.005,0.01),
                        "R" = c(10,20,40))

args <- as.numeric(commandArgs(trailingOnly = TRUE)) # command line arguments
thisrun <- settings[args,] 

N <- sqrt(as.numeric(thisrun[1])) # Number of units
T <- sqrt(as.numeric(thisrun[1]))  # Number of time-periods

noise_sc <- as.numeric(thisrun[2]) # Noise scale 
delta_sc <- 0.1 # delta scale
gamma_sc <- 0.2 # gamma scale
beta_sc <- 0.3 # beta scale
effect_size <- as.numeric(thisrun[3])
R <- as.numeric(thisrun[4])

n.runs <- 1000 # Num. simulation runs

setting <- paste0("N = ", N, ", T = ", T, ", R = ",R, ", noise_sc = ",noise_sc, ", delta_sc = ",delta_sc, ", gamma_sc = ",gamma_sc, ", beta_sc = ", beta_sc, ", effect_size = ", effect_size)
tic(print(paste0("setting: ",setting)))

results <- foreach(i = 1:n.runs, .combine='cbind', .packages =c("MCPanel","matrixStats","boot","Matrix")) %dopar% {
  MCsim(N,T,R,noise_sc,delta_sc,gamma_sc,beta_sc,effect_size,n=i)
}
results <- matrix(unlist(results), ncol = n, byrow = FALSE, .inorder=FALSE) # coerce into matrix
saveRDS(results, paste0("results_","N_",N,"_T_",T,"_R_", R,"_noise_sc_",noise_sc,"_effect_size_",effect_size,"_n_",n,".rds"))

print(toc())

if(doMPI){
  closeCluster(cl) # close down MPIcluster
  mpi.finalize()
}else{
  stopCluster(cl)
}