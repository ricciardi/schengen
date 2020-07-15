######################################################################
# Schengen Simulations #
######################################################################

## Loading Source files
library(MCPanel)
library(glmnet)

# Setup parallel processing 
library(parallel)
library(doParallel)

cores <- parallel::detectCores()
print(paste0('cores registered: ', cores))

cl <- makePSOCKcluster(cores)

doParallel::registerDoParallel(cores) # register cores (<p)

SchengenSim <- function(outcome,sim){
  
  outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  # Use post-treatment (all zeros)
  outcomes.cbw.placebo <- outcomes.cbw
  outcomes.cbw.placebo$mask <- outcomes.cbw$mask[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  outcomes.cbw.placebo$M <- outcomes.cbw$M[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  outcomes.cbw.placebo$W <- outcomes.cbw$W[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  
  Y <- outcomes.cbw.placebo$M # NxT 
  treat <- outcomes.cbw.placebo$mask # NxT masked matrix 
  
  weights <- outcomes.cbw.placebo$W
  weights <- weights[rownames(weights) %in% row.names(Y),]
  weights <- weights[row.names(Y),]  # ensure correct order
  
  weights <- (weights)/(1-weights) 
  
  ## Setting up the configuration
  N <- nrow(treat)
  T <- ncol(treat)
  T0 <- round(c(ncol(outcomes.cbw.placebo$mask)-1, ncol(outcomes.cbw.placebo$mask)/1.25, ncol(outcomes.cbw.placebo$mask)/1.5, ncol(outcomes.cbw.placebo$mask)/2))
  N_t <- ceiling(N*0.5) # no. treated units desired <=N
  num_runs <- 200
  is_simul <- sim ## Whether to simulate Simultaneus Adoption or Staggered Adoption

  ## Matrices for saving RMSE values
  
  MCPanel_RMSE_test <- matrix(0L,num_runs,length(T0))
  DID_RMSE_test <- matrix(0L,num_runs,length(T0))
  ADH_RMSE_test <- matrix(0L,num_runs,length(T0))
  EN_RMSE_test <- matrix(0L,num_runs,length(T0))
  ENT_RMSE_test <- matrix(0L,num_runs,length(T0))
  
  ## Run different methods
  
  for(i in c(1:num_runs)){
    print(paste0(paste0("Run number ", i)," started"))
    ## Fix the treated units in the whole run for a better comparison
    treat_indices <- sort(sample(1:N, N_t))
    for (j in c(1:length(T0))){
      t0 <- T0[j]
      ## Simultaneuous (simul_adapt) or Staggered adoption (stag_adapt)
      if(is_simul == 1){
        treat_mat <- simul_adapt(Y, N_t, t0, treat_indices) 
      }else{
        treat_mat <- stag_adapt(Y, N_t, t0, treat_indices) 
      }
      
      rotate <- function(x) t(apply(x, 2, rev))
      
      treat_mat <- rotate(rotate(treat_mat)) # retrospective analysis
      
      Y_obs <- Y * treat_mat # treated are 0 

      ## -----
      ## ADH
      ## -----
      est_model_ADH <- adh_mp_rows(Y_obs, treat_mat, rel_tol = 0.001)
      est_model_ADH_msk_err <- (est_model_ADH - Y)*(1-treat_mat)
      est_model_ADH_test_RMSE <- sqrt((1/sum(1-treat_mat)) * sum(est_model_ADH_msk_err^2, na.rm = TRUE))
      ADH_RMSE_test[i,j] <- est_model_ADH_test_RMSE
      print(paste("ADH RMSE:", round(est_model_ADH_test_RMSE,3),"run",i))
      

      ## ------
      ## MC-NNM
      ## ------
      
      est_model_MCPanel <- mcnnm_cv(M = Y_obs, mask = treat_mat, W = weights, to_estimate_u = 1, to_estimate_v = 1, num_lam_L = 3, niter = 1000, rel_tol = 1e-03, cv_ratio = 0.8, num_folds = 2, is_quiet = 1) 
      est_model_MCPanel$Mhat <- est_model_MCPanel$L  + replicate(T,est_model_MCPanel$u) + t(replicate(N,est_model_MCPanel$v)) # use X with imputed endogenous values
      est_model_MCPanel$msk_err <- (est_model_MCPanel$Mhat - Y)*(1-treat_mat)
      est_model_MCPanel$test_RMSE <- sqrt((1/sum(1-treat_mat)) * sum(est_model_MCPanel$msk_err^2, na.rm = TRUE))
      MCPanel_RMSE_test[i,j] <- est_model_MCPanel$test_RMSE
      print(paste("MC-NNM RMSE:", round(est_model_MCPanel$test_RMSE,3),"run",i))
      
      ## -----
      ## DID
      ## -----
      
      est_model_DID <- t(DID(t(Y_obs), t(treat_mat)))
      est_model_DID_msk_err <- (est_model_DID - Y)*(1-treat_mat)
      est_model_DID_test_RMSE <- sqrt((1/sum(1-treat_mat)) * sum(est_model_DID_msk_err^2, na.rm = TRUE))
      DID_RMSE_test[i,j] <- est_model_DID_test_RMSE
      print(paste("DID RMSE:", round(est_model_DID_test_RMSE,3),"run",i))
      
      ## -----
      ## VT-EN 
      ## -----
      
      est_model_ENT <- t(en_mp_rows(t(Y_obs), t(treat_mat), num_folds = 3))
      est_model_ENT_msk_err <- (est_model_ENT - Y)*(1-treat_mat)
      est_model_ENT_test_RMSE <- sqrt((1/sum(1-treat_mat)) * sum(est_model_ENT_msk_err^2, na.rm = TRUE))
      ENT_RMSE_test[i,j] <- est_model_ENT_test_RMSE
      print(paste("VT-EN RMSE:", round(est_model_ENT_test_RMSE,3),"run",i))
    }
  }
  
  ## Computing means and standard errors
  
  MCPanel_avg_RMSE <- apply(MCPanel_RMSE_test,2,mean)
  MCPanel_std_error <- apply(MCPanel_RMSE_test,2,sd)/sqrt(num_runs)
  
  DID_avg_RMSE <- apply(DID_RMSE_test,2,mean)
  DID_std_error <- apply(DID_RMSE_test,2,sd)/sqrt(num_runs)
  
  ADH_avg_RMSE <- apply(ADH_RMSE_test,2,mean)
  ADH_std_error <- apply(ADH_RMSE_test,2,sd)/sqrt(num_runs)
  
  ENT_avg_RMSE <- apply(ENT_RMSE_test,2,mean)
  ENT_std_error <- apply(ENT_RMSE_test,2,sd)/sqrt(num_runs)
  
  ## Saving data
  
  df1 <-
    data.frame(
      "y" =  c(DID_avg_RMSE,MCPanel_avg_RMSE,ADH_avg_RMSE,ENT_avg_RMSE),
      "lb" = c(DID_avg_RMSE - 1.96*DID_std_error,
               MCPanel_avg_RMSE - 1.96*MCPanel_std_error, 
               ADH_avg_RMSE - 1.96*ADH_std_error,
               ENT_avg_RMSE - 1.96*ENT_std_error),
      "ub" = c(DID_avg_RMSE + 1.96*DID_std_error, 
               MCPanel_avg_RMSE + 1.96*MCPanel_std_error, 
               ADH_avg_RMSE + 1.96*ADH_std_error,
               ENT_avg_RMSE + 1.96*ENT_std_error),
      "x" = T0/T,
      "Method" = c(replicate(length(T0),"DID"), 
                   replicate(length(T0),"MC-NNM"), 
                   replicate(length(T0),"SCM"),
                   replicate(length(T0),"Vertical")))
  
  filename<-paste0(paste0(paste0(paste0(paste0(paste0(gsub("\\.", "_", o),"_N_", N),"_T_", T),"_numruns_", num_runs), "_num_treated_", N_t), "_simultaneuous_", is_simul),".rds")
  save(df1, file = paste0("results/",filename))

}

# Read data
outcome.vars <- c("CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_3more")

for(o in outcome.vars){
  SchengenSim(outcome=o, sim=1)
}