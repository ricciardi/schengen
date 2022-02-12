#######################################
# placebo simulations on actual data  #
#######################################

## Libraries
library(MCPanel)
library(matrixStats)
library(Matrix)
library(tictoc)
library(data.table)
library(reshape)
library(reshape2)
library(emfactor)
library(parallel)

source('utils.R')
source('IFE.R')

## Setup parallel processing 
# Setup parallel processing

cores <- parallel::detectCores()
print(paste0("number of cores used: ", cores))

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
  library(doParallel)
  library(foreach)
  
  cl <- parallel::makeCluster(cores, outfile="")
  
  doParallel::registerDoParallel(cl) # register cluster
}

SchengenSim <- function(T0,N_t,sim,n,outcome){
  
  # set the seed
  print(paste0("run number: ", n))
  set.seed(n, "L'Ecuyer-CMRG")
  
  # read the data
  outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",outcome,".rds"))
  
  # Use post-treatment (mask is all zeros)
  outcomes.cbw.placebo <- outcomes.cbw
  outcomes.cbw.placebo$mask <- outcomes.cbw$mask[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  outcomes.cbw.placebo$M <- outcomes.cbw$M[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  outcomes.cbw.placebo$W <- outcomes.cbw$W[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  outcomes.cbw.placebo$X <- outcomes.cbw$X[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  outcomes.cbw.placebo$X.hat <- outcomes.cbw$X.hat[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  
  Y <- outcomes.cbw.placebo$M # NxT 
  X <- outcomes.cbw.placebo$X # NxT 
  X.hat <- outcomes.cbw.placebo$X # NxT 
  treat <- outcomes.cbw.placebo$mask # NxT masked matrix 
  
  W <- outcomes.cbw.placebo$W
  W <- W[rownames(W) %in% row.names(Y),]
  W <- W[row.names(Y),]  # ensure correct order

  ## Setting up the configuration
  N <- nrow(treat)
  T <- ncol(treat)
  t0 <- ceiling(T*T0)
  n_t <- ceiling(N*N_t) # no. treated units desired <=N
  is_simul <- sim
  
  ## Run different methods
  
  ## Fix the treated units in the whole run for a better comparison
  treat_indices <- sort(sample(1:N, n_t))
  
  ## Simultaneuous (simul_adapt) or Staggered adoption (stag_adapt)
  if(is_simul == 1){
    treat_mat <- simul_adapt(Y, n_t, t0, treat_indices) 
  }else{
    treat_mat <- stag_adapt(Y, n_t, t0, treat_indices) 
  }
  
  treat_mat <- treat_mat[,c(T:1)] # retrospective analysis
  
  Y_obs <- Y * treat_mat # treated are 0
  
  rownames(Y_obs) <- 1:nrow(Y_obs)
  
  # get vector of initial treatment periods
  
  A <- aggregate(col ~ row,
                 data = which(treat_mat == 1, arr.ind = T),
                 FUN = function(x) x[1])$col
  A[which(A==1)] <- Inf
  
  ST <- which(!is.infinite(A)) # switch treated indices
  AT <- which(is.infinite(A)) # always-treated indices
  
  ## ------ ------ ------ ------ ------
  ## MC-NNM + weights + covariate
  ## ------ ------ ------ ------ ------
  
  est_mc_weights_covars <- mcnnm_wc(M = Y_obs, C = X, mask = treat_mat, W = W, to_estimate_u = 1, to_estimate_v = 1, lambda_L = 0.00108901, lambda_B = 0.100334, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
  while(any(is.na(est_mc_weights_covars$B))){
    est_mc_weights_covars <- mcnnm_wc(M = Y_obs, C = X, mask = treat_mat, W = W, to_estimate_u = 1, to_estimate_v = 1, lambda_L = 0.00108901, lambda_B = 0.100334, niter = 1000, rel_tol = 1e-05, is_quiet = 1)[[1]]
  }
  est_mc_weights_covars$Mhat <- est_mc_weights_covars$L + X.hat*mean(est_mc_weights_covars$B) + replicate(T,est_mc_weights_covars$u) + t(replicate(N,est_mc_weights_covars$v)) # use X with imputed endogenous values
  est_mc_weights_covars$msk_err <- (est_mc_weights_covars$Mhat - Y)*(1-treat_mat)
  est_mc_weights_covars$RMSE <- sqrt((1/sum(1-treat_mat)) * sum(est_mc_weights_covars$msk_err^2, na.rm = TRUE))
  print(paste("MC-NNM (weights +covars) RMSE:", round(est_mc_weights_covars$RMSE,5)))
  
  est_mc_weights_covars$att <- apply(est_mc_weights_covars$msk_err,1,nzmean)[ST]
  est_mc_weights_covars$att.bar <- mean(est_mc_weights_covars$att)
  est_mc_weights_covars$abs.bias <- abs(est_mc_weights_covars$att.bar-0)
  print(paste("MC-NNM (weights + covars) abs. bias:", round(est_mc_weights_covars$abs.bias,5)))
  
  # bootstrap variance estimation
  df_mc_weights_covars <- widetoLong(Y= Y_obs, mask = treat_mat, X = X)
  est_mc_weights_covars$boot_var <- clustered_bootstrap(current_data_realized_long=df_mc_weights_covars, estimator="mc_weights", N=N, T=T, B = 999, est_weights = TRUE, ncores=cores) # bootstrap without covariate
  print(paste("MC-NNM (weights + covars) variance:", round(est_mc_weights_covars$boot_var,5)))
  
  est_mc_weights_covars$cp <- CI_test(est_coefficent=est_mc_weights_covars$att.bar, real_coefficent=0, est_var=est_mc_weights_covars$boot_var)
  print(paste("MC-NNM (weights + covars) CP:", round(est_mc_weights_covars$cp,5)))
  
  ## -----
  ## ADH
  ## -----
  est_model_ADH <- list()
  est_model_ADH$Mhat <- adh_mp_rows(Y_obs, treat_mat)
  est_model_ADH$msk_err <- (est_model_ADH$Mhat - Y)*(1-treat_mat)
  est_model_ADH$RMSE <- sqrt((1/sum(1-treat_mat)) * sum(est_model_ADH$msk_err^2, na.rm = TRUE))
  print(paste("ADH RMSE:", round(est_model_ADH$RMSE,5)))
  
  est_model_ADH$att <- apply(est_model_ADH$msk_err,1,nzmean)[ST]
  est_model_ADH$att.bar <- mean(est_model_ADH$att)
  est_model_ADH$abs.bias <- abs(est_model_ADH$att.bar-0)
  print(paste("ADH abs. bias:", round(est_model_ADH$abs.bias,5)))
  
  # bootstrap variance estimation
  df_compare <- widetoLong(Y= Y_obs, mask = treat_mat, X = NULL)
  est_model_ADH$boot_var <- clustered_bootstrap(current_data_realized_long=df_compare, estimator="ADH", N=N, T=T, B = 999, est_weights = FALSE, ncores=cores)
  print(paste("ADH variance:", round(est_model_ADH$boot_var,5)))
  
  est_model_ADH$cp <- CI_test(est_coefficent=est_model_ADH$att.bar, real_coefficent=0, est_var=est_model_ADH$boot_var)
  print(paste("ADH CP:", round(est_model_ADH$cp,5)))
  
  ## -----
  ## DID
  ## -----
  
  est_model_DID <- list()
  est_model_DID$Mhat<- DID(Y_obs, treat_mat)
  est_model_DID$msk_err <- (est_model_DID$Mhat - Y)*(1-treat_mat)
  est_model_DID$RMSE <- sqrt((1/sum(1-treat_mat)) * sum(est_model_DID$msk_err^2, na.rm = TRUE))
  print(paste("DID RMSE:", round(est_model_DID$RMSE,5)))
  
  est_model_DID$att <- apply(est_model_DID$msk_err,1,nzmean)[ST]
  est_model_DID$att.bar <- mean(est_model_DID$att)
  est_model_DID$abs.bias <- abs(est_model_DID$att.bar-0)
  print(paste("DID abs. bias:", round(est_model_DID$abs.bias,5)))
  
  # bootstrap variance estimation
  est_model_DID$boot_var <- clustered_bootstrap(current_data_realized_long=df_compare, estimator="DID", N=N, T=T, B = 999, est_weights = FALSE, ncores=cores)
  print(paste("DID variance:", round(est_model_DID$boot_var,5)))
  
  est_model_DID$cp <- CI_test(est_coefficent=est_model_DID$att.bar, real_coefficent=0, est_var=est_model_DID$boot_var)
  print(paste("DID CP:", round(est_model_DID$cp,5)))
  
  ## -----
  ## IFE
  ## -----
  
  est_model_IFE <- list()
  est_model_IFE$Mhat <- IFE(Y_obs, treat_mat, k=2)
  est_model_IFE$msk_err <- (est_model_IFE$Mhat - Y)*(1-treat_mat)
  est_model_IFE$RMSE <- sqrt((1/sum(1-treat_mat)) * sum(est_model_IFE$msk_err^2, na.rm = TRUE))
  print(paste("IFE RMSE:", round(est_model_IFE$RMSE,5)))
  
  est_model_IFE$att <- apply(est_model_IFE$msk_err,1,nzmean)[ST]
  est_model_IFE$att.bar <- mean(est_model_IFE$att)
  est_model_IFE$abs.bias <- abs(est_model_IFE$att.bar-0)
  print(paste("IFE abs. bias:", round(est_model_IFE$abs.bias,5)))
  
  # bootstrap variance estimation
  est_model_IFE$boot_var <- clustered_bootstrap(current_data_realized_long=df_compare, estimator="IFE", N=N, T=T, B = 999, est_weights = FALSE, ncores=cores)
  print(paste("IFE variance:", round(est_model_IFE$boot_var,5)))
  
  est_model_IFE$cp <- CI_test(est_coefficent=est_model_IFE$att.bar, real_coefficent=0, est_var=est_model_IFE$boot_var)
  print(paste("IFE CP:", round(est_model_IFE$cp,5)))
  
  gc()
  cat(paste("Done with simulation run number",n, "\n"))
  return(list("N"=N, "T"=T, "T0"=T0, "N_T"=N_t, "outcome"= outcome, "is_simul", is_simul,
              "est_mc_weights_covars_RMSE"=est_mc_weights_covars$RMSE,"est_mc_weights_covars_abs_bias"=est_mc_weights_covars$abs.bias,"est_mc_weights_covars_cp"=est_mc_weights_covars$cp,"est_mc_weights_covars_boot_var"=est_mc_weights_covars$boot_var,
              "est_model_ADH_RMSE"=est_model_ADH$RMSE,"est_model_ADH_abs_bias"=est_model_ADH$abs.bias,"est_model_ADH_cp"=est_model_ADH$cp,"est_model_ADH_boot_var"=est_model_ADH$boot_var,
              "est_model_DID_RMSE"=est_model_DID$RMSE,"est_model_DID_abs_bias"=est_model_DID$abs.bias,"est_model_DID_cp"=est_model_DID$cp,"est_model_DID_boot_var"=est_model_DID$boot_var,
              "est_model_IFE_RMSE"=est_model_IFE$RMSE,"est_model_IFE_abs_bias"=est_model_IFE$abs.bias,"est_model_IFE_cp"=est_model_IFE$cp,"est_model_IFE_boot_var"=est_model_IFE$boot_var))
}

# define settings for simulation
settings <- expand.grid("T0" =c(0.05, 0.25, 0.5),
                        "N_t" =c(0.5))

args <- commandArgs(trailingOnly = TRUE) # command line arguments
thisrun <- settings[as.numeric(args[1]),] 

T0 <- as.numeric(thisrun[1])  # initial time period in staggered adoption setting
N_t <- as.numeric(thisrun[2]) # Number of ST units

is_simul <- 0
outcome <- "CBWbordEMPL"

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

setting <- paste0("T0 = " ,T0, "N_t", N_t)
tic(print(paste0("setting: ",setting)))

results <- foreach(i = 1:n.runs, .combine='rbind', .packages =c("MCPanel","matrixStats","Matrix","data.table","reshape","reshape2","emfactor"), .verbose = FALSE) %dopar% {
  SchengenSim(T0,N_t,sim=is_simul,n=i,outcome="CBWbordEMPL")
}
results
saveRDS(results, paste0(output_dir,"placebo_results_","T0_",T0, "_N_t_", N_t, "_outcome_", outcome,"_is_simul_", is_simul,"_n_",n.runs,".rds"))

print(toc())

if(doMPI){
  closeCluster(cl) # close down MPIcluster
  mpi.finalize()
}else{
  stopCluster(cl)
}