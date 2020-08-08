## Load aggregated data

library(dplyr)
library(readstata13)
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
## N_CBWbord: number of cross-border workers in the region
## CBWbord: share of residents working in another country, which shares the border with the region of residence, unconditional on employment
## CBWbordEMPL: share of residents working in another country, which shares the border with the region of residence, conditional on employment
## average total working hours (Thwusual)

outcomes <- c("CBWbord","CBWbordEMPL","Thwusual")

## Covariates:

covariates <- c("AV_age_22","AV_age_27","AV_age_32","AV_age_37","AV_age_42","AV_age_47","AV_age_52","AV_age_57", # avail for all analyses
                "AV_women",
                "AV_lowEDU","AV_mediumEDU","AV_highEDU",
                "AV_migr",
                "HHincome_COUNTRY",
                "pop")
#covariates.cbw <- c("GDPcapitaR","R_GDPcapitaR","density", "R_HHincome_COUNTRY2","R_HHincome_COUNTRY","XCGvsEURO_05") # analysis specific covars
covariates.cbw <- c("GDPcapitaR","R_GDPcapitaR") # analysis specific covars

# Analyses: 
## ST vs AT (retrospective, X=CBW) 
cbw <- unique(data$REGION[which(data$treated_CBW!=-1)])

always.treated.cbw <- unique(data$REGION[which(data$treated_CBW=="Controls (always-treated)")])
switch.treated.cbw <- unique(data$REGION[which(data$treated_CBW %in% c("2008: Schengen; 2011: FoM", "2007: FoM; 2009: Schengen"))])
length(always.treated.cbw[!always.treated.cbw%in%switch.treated.cbw]) == length(always.treated.cbw) # ensure these groups don't overlap

# Clusters
## For X in (CBW)
## Eastern cluster (treated_X=1)
### SCHENGEN_X=1 from 20081 and FoM_X=1 from 20111
## Swiss cluster (treated_X=2)
### FoM_X=1 from 2007Q2 and SCHENGEN_X=1 from 2009Q1

eastern.cluster.cbw <- unique(data$REGION[which(data$treated_CBW %in% c("2008: Schengen; 2011: FoM"))])

swiss.cluster.cbw <- unique(data$REGION[which(data$treated_CBW %in% c("2007: FoM; 2009: Schengen"))])

# create yearquarter
data$quarter <-NA
data$quarter <- rep(1:4, nrow(data)/4)
data$yearquarter <- as.numeric(paste0(data$year,data$quarter))

# Covariates matrices

covars <- lapply(c(covariates,covariates.cbw), function(c){
  mat <- reshape(data.frame(data[c("REGION","yearquarter",c)]), idvar = "REGION", timevar = "yearquarter", direction = "wide") # N x T
  colnames(mat) <- sub(paste0(c,"."),"", colnames(mat))
  rownames(mat) <- mat$REGION
  mat <- mat[,-1]
})

names(covars) <- c(covariates,covariates.cbw)

# Covariates X_it

covars.cbw <- lapply(c(covariates.cbw), function(i){ 
  subset <- covars[[i]][rownames(covars[[i]]) %in% cbw,]
  return(subset)
}) # N x # predictors

names(covars.cbw ) <- c(covariates.cbw)

for(o in outcomes){
  
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
  
  # Masked matrix for which 0=control units and treated units before treatment and 1=treated units after treatment
  mask.cbw <- matrix(0, nrow = nrow(data.cbw ), 
                             ncol= ncol(data.cbw),
                             dimnames = list(rownames(data.cbw), colnames(data.cbw))) # (N x T)
  
  for(i in switch.treated.cbw){
    mask.cbw[,colnames(mask.cbw)%in%data$yearquarter[data$REGION==i & data$SCHENGEN_CBW==1 & data$FoM_CBW==1]][rownames(mask.cbw)%in%c(i),] <- 1 # retrospective analysis: estimate Y(1)_LT,pre
  }
   mask.cbw[rownames(mask.cbw)%in%rownames(mask.cbw)[rownames(mask.cbw)%in%switch.treated.cbw],] <- abs(mask.cbw[rownames(mask.cbw)%in%rownames(mask.cbw)[rownames(mask.cbw)%in%switch.treated.cbw],]-1) # retrospective analysis: estimate Y(1)_LT,pre

   ## Impute missing & endogenous values of covariates
   
   source('MCEst.R')
   
   if(o=="CBWbordEMPL"){
     best.var.outcome.cbw.m <- as.matrix(covars.cbw[["R_GDPcapitaR"]][rownames(covars.cbw[["R_GDPcapitaR"]])%in%rownames(mask.cbw),])
   } else{
     best.var.outcome.cbw.m <- as.matrix(covars.cbw[["GDPcapitaR"]][rownames(covars.cbw[["GDPcapitaR"]])%in%rownames(mask.cbw),])
   }
   
   # Masked matrix for which 0=observed, 1 missing
   mask.cbw.missing <- matrix(0, nrow = nrow(data.cbw ), 
                              ncol= ncol(data.cbw),
                              dimnames = list(rownames(data.cbw), colnames(data.cbw))) # (N x T)
   mask.cbw.missing[which(is.na(best.var.outcome.cbw.m))] <- 1
   best.var.outcome.cbw.m[is.na(best.var.outcome.cbw.m)] <- 0 # NA's are 0 in outcome matrix
   
   outcomes.impute.cbw <- list("M"=best.var.outcome.cbw.m, 
                               "mask"=mask.cbw.missing, 
                               "W"= matrix(0.5, nrow(mask.cbw),ncol(mask.cbw),
                                           dimnames = list(rownames(mask.cbw), colnames(mask.cbw))))# weights are equal
   
   impute.best.var.outcome.cbw <- MCEst(outcomes.impute.cbw, rev=TRUE, covars=FALSE, fe=FALSE)
   best.var.outcome.cbw.m.imputed <- best.var.outcome.cbw.m*(1-mask.cbw.missing) + impute.best.var.outcome.cbw$Mhat*mask.cbw.missing # only missing values imputed
   
   outcomes.impute.endog.cbw <- list("M"=best.var.outcome.cbw.m.imputed, # imputed data
                               "mask"=mask.cbw, 
                               "W"= matrix(0.5, nrow(mask.cbw),ncol(mask.cbw),
                                           dimnames = list(rownames(mask.cbw), colnames(mask.cbw))))# weights are equal
   
   impute.endog.best.var.outcome.cbw <- MCEst(outcomes.impute.endog.cbw, rev=TRUE, covars=FALSE, fe=FALSE)
   best.var.outcome.cbw.hat <- best.var.outcome.cbw.m.imputed*(1-mask.cbw) + impute.endog.best.var.outcome.cbw$Mhat*mask.cbw  # only endogenous values imputed
   
   colnames(best.var.outcome.cbw.hat) <- colnames(mask.cbw)
   rownames(best.var.outcome.cbw.hat) <- rownames(mask.cbw)
   
  ## Estimate propensity scores by matrix completion
   
   propensity.model.cbw.data <- list("M"=mask.cbw, 
                                     "X" = best.var.outcome.cbw.m.imputed, # var with endogenous values + imputed missing values
                                     "X.hat"= best.var.outcome.cbw.hat, # var with imputed endogenous values
                                     "mask" = matrix(0, nrow(mask.cbw),ncol(mask.cbw),
                                                     dimnames = list(rownames(mask.cbw), colnames(mask.cbw))), # no missing entries
                                     "W"= matrix(0.5, nrow(mask.cbw),ncol(mask.cbw),
                                                 dimnames = list(rownames(mask.cbw), colnames(mask.cbw)))) 
   
   propensity.model.cbw <- MCEst(propensity.model.cbw.data, rev=TRUE, covars=TRUE, fe=TRUE) 
   
   propensity.model.cbw.values <- propensity.model.cbw$Mhat # L + X_hat*B + u + v
   
   colnames(propensity.model.cbw.values) <- colnames(mask.cbw)
   rownames(propensity.model.cbw.values) <- rownames(mask.cbw)
   
  ## Elapsed time Weights (future and past)
  
  z.cbw.eastern <- round(c(seq(0.999, 0.7, length.out=which(colnames(mask.cbw)=="20111")),
                           seq(0.712, 0.999, length.out=ncol(mask.cbw)-which(colnames(mask.cbw)=="20111"))),3) # "earliest" combined treatment
  
  z.cbw.swiss <- round(c(seq(0.999, 0.7, length.out=which(colnames(mask.cbw)=="20091")),
                         seq(0.719, 0.999, length.out=ncol(mask.cbw)-which(colnames(mask.cbw)=="20091"))),3)
  
  p.weights.cbw <- matrix(propensity.model.cbw.values%*%diag(z.cbw.eastern), 
                          nrow = nrow(data.cbw ), 
                          ncol= ncol(data.cbw),
                          dimnames = list(rownames(data.cbw), colnames(data.cbw)), byrow = TRUE) # (N x T)
  
  p.weights.cbw[rownames(p.weights.cbw) %in% swiss.cluster.cbw,] <- matrix(propensity.model.cbw.values[rownames(propensity.model.cbw.values) %in% swiss.cluster.cbw,]%*%diag(z.cbw.swiss), 
                                                                           nrow = length(swiss.cluster.cbw), 
                                                                           ncol= ncol(data.cbw), byrow = TRUE) # (N x T)
  # Save
  
  outcomes.cbw <- list("M"=data.cbw, "mask"=mask.cbw, "W"= p.weights.cbw,"W.no.z"= propensity.model.cbw.values,
                       "X"=best.var.outcome.cbw.m.imputed, "X.hat"=best.var.outcome.cbw.hat,
                       "mc.outcome"=impute.best.var.outcome.cbw, "mc.propensity"=propensity.model.cbw,
                               "treated"=rownames(mask.cbw)[rownames(mask.cbw)%in%switch.treated.cbw],
                               "control"=rownames(mask.cbw)[rownames(mask.cbw)%in%always.treated.cbw],
                               "eastern"=eastern.cluster.cbw, "swiss"=swiss.cluster.cbw)
  
  saveRDS(outcomes.cbw, paste0("data/outcomes-cbw-",o,".rds"))
}