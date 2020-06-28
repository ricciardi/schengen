## Load aggregated data

library(dplyr)
library(readstata13)
library(glmnet)
library(caret)
library(MCPanel)

# Setup parallel processing 
library(parallel)
library(doParallel)

cores <- detectCores()

cl <- parallel::makeForkCluster(cores)

doParallel::registerDoParallel(cores) # register cores (<p)

RNGkind("L'Ecuyer-CMRG") # ensure random number generation

data <- read.dta13("FINAL.dta", generate.factors=T) # includes variables for 2 analyses + covariates 

## Outcomes:
## CBWbord: share of residents working in another country, which shares the border with the region of residence, unconditional on employment
## CBWbordEMPL: share of residents working in another country, which shares the border with the region of residence, conditional on employment
## regional employment rate (empl),
## average total working hours (Thwusual)
## unemployed rate (unempl)
## inactivity rate (inact)
## % of unemployed with unemployment duration less than 1 month (seekdur_0), 1-2 months (seekdur_1_2), 3 months or more (seekdur_3more).

outcomes <- c("CBWbord","CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_0","seekdur_1_2","seekdur_3more")

## Covariates:

covariates <- c("AV_age_22","AV_age_27","AV_age_32","AV_age_37","AV_age_42","AV_age_47","AV_age_52","AV_age_57", # avail for all analyses
                "AV_women",
                "AV_lowEDU","AV_mediumEDU","AV_highEDU",
                "AV_migr",
                "HHincome_COUNTRY")
covariates.cbw <- c("GDPcapitaR","R_GDPcapitaR","density", "R_HHincome_COUNTRY2","R_HHincome_COUNTRY","XCGvsEURO_05") # analysis specific covars
covariates.lm <- c("AV_HHsize","AV_HHsizeCHILD","AV_HHsizeOLD","AV_single","AV_oneadultNOCHILD","AV_oneadultCHILD",   # drop Lang b.c it is time invariant
                "AV_coupleNOCHILD","AV_coupleCHILD","AV_pcWORK")

# Analyses: 
## ST vs AT (retrospective, X=CBW) 
cbw <- unique(data$REGION[which(data$treated_CBW!=-1)])

always.treated.cbw <- unique(data$REGION[which(data$treated_CBW=="Controls (always-treated)")])
switch.treated.cbw <- unique(data$REGION[which(data$treated_CBW %in% c("2008: Schengen; 2011: FoM", "2007: FoM; 2009: Schengen"))])
length(always.treated.cbw[!always.treated.cbw%in%switch.treated.cbw]) == length(always.treated.cbw) # ensure these groups don't overlap

## ST vs NT (forward, X=LM) 
lm <- unique(data$REGION[which(data$treated_LM!=-1)])

never.treated.lm <- unique(data$REGION[which(data$treated_LM=="Controls (internal regions)")])
switch.treated.lm <- unique(data$REGION[which(data$treated_LM==c("2008: Schengen; 2011: FoM", "2007: FoM; 2009: Schengen"))]) # Slovenia (country=="SI") are only in the CBW analysis since they only have treated regions and no suitable internal control regions for the LM analysis
length(never.treated.lm[!never.treated.lm%in%switch.treated.lm]) == length(never.treated.lm) # ensure these groups don't overlap

# Clusters
## For X in (CBW,LM)
## Eastern cluster (treated_X=1)
### SCHENGEN_X=1 from 20081 and FoM_X=1 from 20111
## Swiss cluster (treated_X=2)
### FoM_X=1 from 2007Q2 and SCHENGEN_X=1 from 2009Q1

eastern.cluster.cbw <- unique(data$REGION[which(data$treated_CBW %in% c("2008: Schengen; 2011: FoM"))])
eastern.cluster.lm <- unique(data$REGION[which(data$treated_LM %in% c("2008: Schengen; 2011: FoM"))])

swiss.cluster.cbw <- unique(data$REGION[which(data$treated_CBW %in% c("2007: FoM; 2009: Schengen"))])
swiss.cluster.lm <- unique(data$REGION[which(data$treated_LM %in% c("2007: FoM; 2009: Schengen"))])

# create yearquarter
data$quarter <-NA
data$quarter <- rep(1:4, nrow(data)/4)
data$yearquarter <- as.numeric(paste0(data$year,data$quarter))

# Covariates matrices

covars <- lapply(c(covariates,covariates.cbw, covariates.lm), function(c){
  mat <- reshape(data.frame(data[c("REGION","yearquarter",c)]), idvar = "REGION", timevar = "yearquarter", direction = "wide") # N x T
  colnames(mat) <- sub(paste0(c,"."),"", colnames(mat))
  rownames(mat) <- mat$REGION
  mat <- mat[,-1]
})

names(covars) <- c(covariates,covariates.cbw, covariates.lm)

# Covariates X_it

covars.cbw <- lapply(c(covariates,covariates.cbw), function(i){ 
  subset <- covars[[i]][rownames(covars[[i]]) %in% cbw,]
  impute <- predict(preProcess(subset, method = c("medianImpute")), subset) # Impute missing with column medians
  return(impute)
}) # N x # predictors

names(covars.cbw ) <- c(covariates,covariates.cbw)

covars.lm <- lapply(c(covariates,covariates.lm), function(i){ 
  subset <- covars[[i]][rownames(covars[[i]]) %in% lm,]
  impute <- predict(preProcess(subset, method = c("medianImpute")), subset) # Impute missing with column medians
  return(impute)
}) # N x # predictors

names(covars.lm) <- c(covariates,covariates.lm)

for(o in outcomes){
  print(o)
  
  # Outcomes matrix
  data.m <- reshape(data.frame(data[c("REGION","yearquarter",o)]), idvar = "REGION", timevar = "yearquarter", direction = "wide") # N x T
  colnames(data.m) <- sub(paste0(o,"."),"", colnames(data.m))
  rownames(data.m) <- data.m$REGION
  data.m <- data.m[,-1]
  
  DropVariance <- function(mat){ # Remove units with no variance
    drop <- names(which(apply(t(mat), 2, var) == 0))
    return(as.matrix(mat[!rownames(mat)%in%drop,]))
  }
  
  data.cbw <-  as.matrix(data.m[rownames(data.m) %in% cbw,])
  if(length(names(which(apply(t(data.cbw), 2, var) == 0)))>0){
    data.cbw <- DropVariance(data.cbw)
  } 
  
  data.lm <-  as.matrix(data.m[rownames(data.m) %in% lm,])
  if(length(names(which(apply(t(data.lm), 2, var) == 0)))>0){
    data.lm <- DropVariance(data.lm)
  } 
  
  # Masked matrix for which 0=control units and treated units before treatment and 1=treated units after treatment
  mask.cbw <- matrix(0, nrow = nrow(data.cbw ), 
                             ncol= ncol(data.cbw),
                             dimnames = list(rownames(data.cbw), colnames(data.cbw))) # (N x T)
  
  mask.lm <- matrix(0, nrow = nrow(data.lm ), 
                            ncol= ncol(data.lm),
                            dimnames = list(rownames(data.lm), colnames(data.lm))) # (N x T)
  
  # for the eastern group, the treatment period is 20111-20181 (both FoM and Schengen in place)
  # for the Swiss group, we have 20091-20181 as treatment period (both FoM and Schengen in place)
  for(i in switch.treated.cbw){
    mask.cbw[,colnames(mask.cbw)%in%data$yearquarter[data$REGION==i & data$SCHENGEN_CBW==1 & data$FoM_CBW==1]][rownames(mask.cbw)%in%c(i),] <- 1 # retrospective analysis: estimate Y(1)_LT,pre
  }
   mask.cbw[rownames(mask.cbw)%in%rownames(mask.cbw)[rownames(mask.cbw)%in%switch.treated.cbw],] <- abs(mask.cbw[rownames(mask.cbw)%in%rownames(mask.cbw)[rownames(mask.cbw)%in%switch.treated.cbw],]-1) # retrospective analysis: estimate Y(1)_LT,pre

  for(i in switch.treated.lm){
    mask.lm[,colnames(mask.lm)%in%data$yearquarter[data$REGION==i & data$SCHENGEN_LM==1 & data$FoM_LM==1]][rownames(mask.lm)%in%c(i),] <- 1 
  }
  
  ## “double-lasso” covariate selection procedure (Belloni, et al., 2014)
  
  # CBW
   
  covars.cbw.reduced <- lapply(1:length(covars.cbw), function(i){ # make sure same dimension as mask
    samedims <- covars.cbw[[i]][rownames(covars.cbw[[i]])%in%rownames(mask.cbw),]
    select <- samedims[,which(colnames(samedims)=='20111'):ncol(samedims)] # only posttreatment years
    return(select)
  })
  names(covars.cbw.reduced) <- c(covariates,covariates.cbw)
  
  covars.cbw.combined <- as.matrix(bind_cols(lapply(1:length(covars.cbw.reduced), function(i){ 
    combined <- covars.cbw.reduced[[i]][,!duplicated(t(covars.cbw.reduced[[i]]))] # combine and remove duplicates (quarters)
    colnames(combined) <- paste0(names(covars.cbw.reduced)[i],colnames(combined)) # naming for variable selection
    return(combined)
  }))) # N x # predictors
  rownames(covars.cbw.combined) <- rownames(mask.cbw)

  # Step 1: Fit a lasso regression predicting the dependent variable
   
  cvfit.outcome.cbw <- cv.glmnet(x=covars.cbw.combined, # post-treatment series
                             y=data.cbw[,1:which(colnames(data.cbw)=='20091')], # pre treatment series 
                             family="mgaussian",
                             standardize.response = TRUE,
                             parallel = TRUE)
  
  # variables with non-zero estimated coefficients: 
  tmp_coeffs <- coef(cvfit.outcome.cbw, s = "lambda.min")[[1]] # same nonzero variables for time series
  tmp_coeffs <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)[-1,] # rm intercept
  op <- options(warn=2)
  best.var.outcome.cbw <- try(as.character(tmp_coeffs$name[which(abs(tmp_coeffs$coefficient) == max(abs(tmp_coeffs$coefficient)))])) # select highest nonzero var
  warn.var.outcome.cbw <- FALSE
  
  if(is(best.var.outcome.cbw ,"try-error") || is.null(best.var.outcome.cbw)){ # if all nonzero randomly select covar
    warn.var.outcome.cbw <- TRUE
    best.var.outcome.cbw <- sample(colnames(covars.cbw.combined),1)
  }
  options(op)
  
  # Step 2: Fit a lasso logistic regression predicting treatment

  cvfit.treatment.cbw  <- cv.glmnet(x=covars.cbw.combined,# post-treatment series
                            y=as.factor(mask.cbw[,"20084"]), # 1 if switch.treated.cbw
                            family="binomial", parallel = TRUE)
  
  # variables with non-zero estimated coefficients: 
  tmp_coeffs <- coef(cvfit.treatment.cbw, s = "lambda.min")
  tmp_coeffs <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)[-1,] # rm intercept
  op <- options(warn=2)
  best.var.treatment.cbw <- try(as.character(tmp_coeffs$name[which(abs(tmp_coeffs$coefficient) == max(abs(tmp_coeffs$coefficient)))])) # select highest nonzero var
  warn.var.treatment.cbw <- FALSE
  
  if(is(best.var.treatment.cbw ,"try-error") || is.null(best.var.treatment.cbw)){ # if all nonzero randomly select covar
    warn.var.treatment.cbw <- TRUE
    best.var.treatment.cbw <- sample(colnames(covars.cbw.combined),1)
  }
  options(op)
  
  ## Impute endogenous values of best covars
  
  # Outcome model variable
  source('MCEst.R')
  
  best.var.outcome.cbw.m <- as.matrix(covars.cbw[[gsub('.{5}$', '', best.var.outcome.cbw)]][rownames(covars.cbw[[1]])%in%rownames(mask.cbw),])
  outcomes.impute.cbw <- list("M"=best.var.outcome.cbw.m, # outcome is best covariate for outcome model
                                "mask"=mask.cbw, 
                              "W"= matrix(1-.Machine
                                          $double.eps, nrow(mask.cbw),ncol(mask.cbw),
                                          dimnames = list(rownames(mask.cbw), colnames(mask.cbw))))# weights are virtually 1
  
  impute.best.var.outcome.cbw <- MCEst(outcomes.impute.cbw, rev=TRUE, covars=FALSE, nofes=TRUE) # run with no FEs
  best.var.outcome.cbw.hat <- best.var.outcome.cbw.m*(1-mask.cbw) + impute.best.var.outcome.cbw$Mhat*mask.cbw # only endogenous values imputed

  colnames(best.var.outcome.cbw.hat) <- colnames(mask.cbw)
  rownames(best.var.outcome.cbw.hat) <- rownames(mask.cbw)
  
  # Propensity model variable
  
  best.var.treatment.cbw.m <- as.matrix(covars.cbw[[gsub('.{5}$', '', best.var.treatment.cbw)]][rownames(covars.cbw[[1]])%in%rownames(mask.cbw),])
  treatment.impute.cbw <- list("M"=best.var.treatment.cbw.m, # outcome is best covariate for propensity model
                              "mask"=mask.cbw, 
                              "W"= matrix(1-.Machine
                                          $double.eps, nrow(mask.cbw),ncol(mask.cbw),
                                          dimnames = list(rownames(mask.cbw), colnames(mask.cbw))))# weights are virtually 1
  
  impute.best.var.treatment.cbw <- MCEst(treatment.impute.cbw, rev=TRUE, covars=FALSE, nofes=TRUE) # run with no FEs
  best.var.treatment.cbw.hat <- best.var.treatment.cbw.m*(1-mask.cbw) + impute.best.var.treatment.cbw$Mhat*mask.cbw # only endogenous values imputed
  
  colnames(best.var.treatment.cbw.hat) <- colnames(mask.cbw)
  rownames(best.var.treatment.cbw.hat) <- rownames(mask.cbw)
  
  ## Est propensity scores as a function of best propensity model variable
  
  propensity.model.cbw <- mcnnm_wc_cv(M = mask.cbw, 
                                      C = best.var.treatment.cbw.m, # best propensity model var with endogenous values
                                      mask = matrix(1, nrow(mask.cbw),ncol(mask.cbw),
                                             dimnames = list(rownames(mask.cbw), colnames(mask.cbw))), # no missing entries
                                      W = matrix(1,nrow(mask.cbw),ncol(mask.cbw)),
                                      to_normalize = 1, 
                                      to_estimate_u = 0, to_estimate_v = 0, # no covariates
                                      num_lam_L = 5, num_lam_B = 5, niter = 1000, rel_tol = 1e-03, cv_ratio = 0.8, num_folds = 2, is_quiet = 1)
  
  propensity.model.cbw$E <-propensity.model.cbw$L +  replicate(ncol(mask.cbw), as.vector(best.var.treatment.cbw.hat%*%propensity.model.cbw$B)) # reconstruct with imputed endogenous values
  
  colnames(propensity.model.cbw$E) <- colnames(mask.cbw)
  rownames(propensity.model.cbw$E) <- rownames(mask.cbw)
  
  ## Elapsed time Weights
  
  z.cbw <- round(c(seq(1, 0.7, length.out=which(colnames(mask.cbw)=="20091")),
                           seq(0.719, 1, length.out=ncol(mask.cbw)-which(colnames(mask.cbw)=="20091"))),3) # earliest combined treatment
  
  p.weights.cbw <- matrix(0, nrow = nrow(data.cbw ), 
                                        ncol= ncol(data.cbw),
                                        dimnames = list(rownames(data.cbw), colnames(data.cbw))) # (N x T)
  
  p.weights.cbw <- propensity.model.cbw$E%*%diag(z.cbw)

  # LM
  
  covars.lm.reduced <- lapply(1:length(covars.lm), function(i){ # make sure same dimension as mask
    samedims <- covars.lm[[i]][rownames(covars.lm[[i]])%in%rownames(mask.lm),]
    select <- samedims[,1:which(colnames(samedims)=='20071')] # only pretreatment years
    return(select)
  })
  names(covars.lm.reduced) <- c(covariates,covariates.lm)
  
  covars.lm.combined <- as.matrix(bind_cols(lapply(1:length(covars.lm.reduced), function(i){ 
    combined <- covars.lm.reduced[[i]][,!duplicated(t(covars.lm.reduced[[i]]))] # combine and remove duplicates (quarters)
    colnames(combined) <- paste0(names(covars.lm.reduced)[i],colnames(combined)) # naming for variable selection
    return(combined)
  }))) # N x # predictors
  rownames(covars.lm.combined) <- rownames(mask.lm)

  # Step 1: Fit a lasso regression predicting the dependent variable
  
  cvfit.outcome.lm <- cv.glmnet(x=covars.lm.combined, # pre-treatment series
                                 y=data.lm[,which(colnames(data.lm)=='20111'):ncol(data.lm)], # post-treatment series 
                                 family="mgaussian",
                                 standardize.response = TRUE,
                                 parallel = TRUE)
  
  # variables with non-zero estimated coefficients: 
  tmp_coeffs <- coef(cvfit.outcome.lm, s = "lambda.min")[[1]] # same nonzero variables for time series
  tmp_coeffs <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)[-1,] # rm intercept
  op <- options(warn=2)
  best.var.outcome.lm <- try(as.character(tmp_coeffs$name[which(abs(tmp_coeffs$coefficient) == max(abs(tmp_coeffs$coefficient)))])) # select highest nonzero var
  warn.var.outcome.lm <- FALSE
  
  if(is(best.var.outcome.lm ,"try-error") || is.null(best.var.outcome.lm)){ # if all nonzero randomly select covar
    warn.var.outcome.lm <- TRUE
    best.var.outcome.lm <- sample(colnames(covars.lm.combined),1)
  }
  options(op)
  
  # Step 2: Fit a lasso logistic regression predicting treatment
  
  cvfit.treatment.lm  <- cv.glmnet(x=covars.lm.combined, # pre-treatment series
                                    y=as.factor(mask.lm[,"20111"]), # 1 if switch.treated.lm
                                    family="binomial", parallel = TRUE)
  
  # variables with non-zero estimated coefficients: 
  tmp_coeffs <- coef(cvfit.treatment.lm, s = "lambda.min")
  tmp_coeffs <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)[-1,] # rm intercept
  op <- options(warn=2)
  best.var.treatment.lm <- try(as.character(tmp_coeffs$name[which(abs(tmp_coeffs$coefficient) == max(abs(tmp_coeffs$coefficient)))])) # select highest nonzero var
  warn.var.treatment.lm <- FALSE
  if(is(best.var.treatment.lm ,"try-error") || is.null(best.var.treatment.lm)){ # if all nonzero
    warn.var.treatment.lm <- TRUE
    best.var.treatment.lm <- sample(colnames(covars.lm.combined),1)
  }
  options(op)
  
  ## Impute endogenous values of best covars
  
  # Outcome model variable
  best.var.outcome.lm.m <- as.matrix(covars.lm[[gsub('.{5}$', '', best.var.outcome.lm)]][rownames(covars.lm[[1]])%in%rownames(mask.lm),])
  outcomes.impute.lm <- list("M"=best.var.outcome.lm.m, # outcome is best covariate
                              "mask"=mask.lm, 
                              "W"= matrix(1-.Machine
                                          $double.eps, nrow(mask.lm),ncol(mask.lm),
                                          dimnames = list(rownames(mask.lm), colnames(mask.lm))))# weights are virtually 1
  
  impute.best.var.outcome.lm <- MCEst(outcomes.impute.lm, rev=FALSE, covars=FALSE, nofes=TRUE) # run with no FEs
  best.var.outcome.lm.hat <- best.var.outcome.lm.m*(1-mask.lm) + impute.best.var.outcome.lm$Mhat*mask.lm # only endogenous values imputed
  
  colnames(best.var.outcome.lm.hat) <- colnames(mask.lm)
  rownames(best.var.outcome.lm.hat) <- rownames(mask.lm)
  
  # Propensity model variable
  
  best.var.treatment.lm.m <- as.matrix(covars.lm[[gsub('.{5}$', '', best.var.treatment.lm)]][rownames(covars.lm[[1]])%in%rownames(mask.lm),])
  treatment.impute.lm <- list("M"=best.var.treatment.lm.m, # outcome is best covariate
                               "mask"=mask.lm, 
                               "W"= matrix(1-.Machine
                                           $double.eps, nrow(mask.lm),ncol(mask.lm),
                                           dimnames = list(rownames(mask.lm), colnames(mask.lm))))# weights are virtually 1
  
  impute.best.var.treatment.lm <- MCEst(treatment.impute.lm, rev=FALSE, covars=FALSE, nofes=TRUE) # run with no FEs
  best.var.treatment.lm.hat <- best.var.treatment.lm.m*(1-mask.lm) + impute.best.var.treatment.lm$Mhat*mask.lm # only endogenous values imputed
  
  colnames(best.var.treatment.lm.hat) <- colnames(mask.lm)
  rownames(best.var.treatment.lm.hat) <- rownames(mask.lm)
  
  ## Est propensity scores as a function of best propensity model variable (with endogenous values)
  
  propensity.model.lm <- mcnnm_wc_cv(M = mask.lm, 
                                      C = best.var.treatment.lm.m, # var with endogenous values
                                      mask = matrix(1, nrow(mask.lm),ncol(mask.lm),
                                                    dimnames = list(rownames(mask.lm), colnames(mask.lm))), # no missing entries
                                      W = matrix(1,nrow(mask.lm),ncol(mask.lm)),
                                      to_normalize = 1, 
                                      to_estimate_u = 0, to_estimate_v = 0, 
                                      num_lam_L = 5, num_lam_B = 5, niter = 1000, rel_tol = 1e-03, cv_ratio = 0.8, num_folds = 2, is_quiet = 1)
  
  propensity.model.lm$E <-propensity.model.lm$L +  replicate(ncol(mask.lm), as.vector(best.var.treatment.lm.hat%*%propensity.model.lm$B)) # reconstruct with imputed endogenous values
  
  colnames(propensity.model.lm$E) <- colnames(mask.lm)
  rownames(propensity.model.lm$E) <- rownames(mask.lm)
  
  ## Elapsed time Weights
  
  z.lm <- round(c(seq(1, 0.7, length.out=which(colnames(mask.lm)=="20111")),
                   seq(0.719, 1, length.out=ncol(mask.lm)-which(colnames(mask.lm)=="20111"))),3) # latest combined treatment
  
  p.weights.lm <- matrix(0, nrow = nrow(data.lm ), 
                          ncol= ncol(data.lm),
                          dimnames = list(rownames(data.lm), colnames(data.lm))) # (N x T)
  
  p.weights.lm <- propensity.model.lm$E%*%diag(z.lm)
  
 
  # Save
  
  outcomes.cbw <- list("M"=data.cbw, "mask"=mask.cbw, "W"= p.weights.cbw, "X"=best.var.outcome.cbw.m, "X.hat"=best.var.outcome.cbw.hat,
                       "outcome.X" = best.var.outcome.cbw, "propensity.X" = best.var.treatment.cbw,
                       "outcome.X.warn"=warn.var.outcome.cbw, "propensity.X.warn"=warn.var.treatment.cbw, 
                       "mc.outcome"=impute.best.var.outcome.cbw, "mc.propensity"=propensity.model.cbw, 
                       "z"=z.cbw, "lasso.outcome"=cvfit.outcome.cbw, "lasso.propensity"=cvfit.treatment.cbw,
                               "treated"=rownames(mask.cbw)[rownames(mask.cbw)%in%switch.treated.cbw],
                               "control"=rownames(mask.cbw)[rownames(mask.cbw)%in%always.treated.cbw],
                               "eastern"=eastern.cluster.cbw, "swiss"=swiss.cluster.cbw)
  outcomes.lm <- list("M"=data.lm, "mask"=mask.lm, "W"= p.weights.lm, "X"=best.var.outcome.lm.m, "X.hat"=best.var.outcome.lm.hat,
                      "outcome.X" = best.var.outcome.lm, "propensity.X" = best.var.treatment.lm,
                      "outcome.X.warn"=warn.var.outcome.lm, "propensity.X.warn"=warn.var.treatment.lm, 
                      "mc.outcome"=impute.best.var.outcome.lm, "mc.propensity"=propensity.model.lm, 
                      "z"=z.lm, "lasso.outcome"=cvfit.outcome.lm, "lasso.propensity"=cvfit.treatment.lm,
                              "treated"=rownames(mask.lm)[rownames(mask.lm)%in%switch.treated.lm],
                              "control"=rownames(mask.lm)[rownames(mask.lm)%in%never.treated.lm],
                              "eastern"=eastern.cluster.lm, "swiss"=swiss.cluster.lm)
  
  saveRDS(outcomes.cbw, paste0("data/outcomes-cbw-",o,".rds"))
  saveRDS(outcomes.lm, paste0("data/outcomes-lm-",o,".rds"))
}