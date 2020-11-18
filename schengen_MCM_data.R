## Load aggregated data

library(dplyr)
library(readstata13)
library(caret)
library(MCPanel)
library(Matrix)

# Setup parallel processing 
library(parallel)
library(doParallel)

cores <- detectCores()

cl <- parallel::makeForkCluster(cores)

doParallel::registerDoParallel(cores) # register cores (<p)

RNGkind("L'Ecuyer-CMRG") # ensure random number generation

data <- read.dta13("Collapsed_05_19.dta", generate.factors=T)

## Outcomes:

outcomes <- c("CBWbord","CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_0","seekdur_1_2","seekdur_3more")

## Covariates:

covariates <- c("AV_age_22","AV_age_27","AV_age_32","AV_age_37","AV_age_42","AV_age_47","AV_age_52","AV_age_57", # avail for all analyses
                "AV_women",
                "AV_lowEDU","AV_mediumEDU","AV_highEDU",
                "AV_migr",
                "HHincome_COUNTRY",
                "pop")
covariates.cbw <- c("GDPcapitaR","R_GDPcapitaR") # analysis specific covars

# Analyses: 
## ST vs AT (retrospective, X=CBW) 
cbw <- unique(data$REGION[which(data$treated_CBW!=-1)])

always.treated.cbw <- unique(data$REGION[which(data$treated_CBW==0)]) # "Controls (always-treated)"
switch.treated.cbw <- unique(data$REGION[which(data$treated_CBW %in% c(1:2))]) # "2008: Schengen; 2011: FoM", "2007: FoM; 2009: Schengen"
length(always.treated.cbw[!always.treated.cbw%in%switch.treated.cbw]) == length(always.treated.cbw) # ensure these groups don't overlap

# Clusters
## For X in (CBW)
## Eastern cluster (treated_X=1)
### SCHENGEN_X=1 from 20081 and FoM_X=1 from 20111
## Swiss cluster (treated_X=2)
### FoM_X=1 from 2007Q2 and SCHENGEN_X=1 from 2009Q1

eastern.cluster.cbw <- unique(data$REGION[which(data$treated_CBW==1)]) # "2008: Schengen; 2011: FoM"

swiss.cluster.cbw <- unique(data$REGION[which(data$treated_CBW==2)]) #"2007: FoM; 2009: Schengen"

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

   ## Elapsed time Weights (future and past)
   
   z.cbw.eastern <- c(rev(SSlogis(1:which(colnames(mask.cbw)=="20111"), Asym = 1, xmid = 0.85, scal = 8)),
                      SSlogis(1:(ncol(mask.cbw)-which(colnames(mask.cbw)=="20111")), Asym = 1, xmid = 0.85, scal = 8))
   
   z.cbw.swiss <- c(rev(SSlogis(1:which(colnames(mask.cbw)=="20091"), Asym = 1, xmid = 0.85, scal = 8)),
                    SSlogis(1:(ncol(mask.cbw)-which(colnames(mask.cbw)=="20091")), Asym = 1, xmid = 0.85, scal = 8))
   
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
   
   p.weights.cbw.equal <- matrix(0.5, nrow(mask.cbw),ncol(mask.cbw),
                                 dimnames = list(rownames(mask.cbw), colnames(mask.cbw)))
   
   outcomes.impute.cbw <- list("M"=best.var.outcome.cbw.m, 
                               "mask"=mask.cbw.missing, 
                               "W"= p.weights.cbw.equal, # weights are equal 
                               "z.cbw.eastern"=z.cbw.eastern,
                               "z.cbw.swiss"= z.cbw.swiss,
                               "eastern"=eastern.cluster.cbw, 
                               "swiss"=swiss.cluster.cbw)
   
   impute.best.var.outcome.cbw.eastern <- MCEst(outcomes.impute.cbw, cluster='eastern', rev=TRUE, covars=FALSE)
   impute.best.var.outcome.cbw.swiss <- MCEst(outcomes.impute.cbw, cluster='swiss', rev=TRUE, covars=FALSE)
   best.var.outcome.cbw.m.imputed <- matrix(NA, nrow(mask.cbw),ncol(mask.cbw),
                                            dimnames = list(rownames(mask.cbw), colnames(mask.cbw)))
   best.var.outcome.cbw.m.imputed[!rownames(best.var.outcome.cbw.m.imputed)%in%swiss.cluster.cbw,] <- best.var.outcome.cbw.m[!rownames(best.var.outcome.cbw.m)%in%swiss.cluster.cbw,]*(1-mask.cbw.missing[!rownames(mask.cbw.missing)%in%swiss.cluster.cbw,]) + 
     impute.best.var.outcome.cbw.eastern$Mhat*mask.cbw.missing[!rownames(mask.cbw.missing)%in%swiss.cluster.cbw,] # only missing values imputed 
   
   best.var.outcome.cbw.m.imputed[rownames(best.var.outcome.cbw.m.imputed)%in%swiss.cluster.cbw,] <- (best.var.outcome.cbw.m[!rownames(best.var.outcome.cbw.m)%in%eastern.cluster.cbw,]*(1-mask.cbw.missing[!rownames(mask.cbw.missing)%in%eastern.cluster.cbw,]) + 
     impute.best.var.outcome.cbw.swiss$Mhat*mask.cbw.missing[!rownames(mask.cbw.missing)%in%eastern.cluster.cbw,])[swiss.cluster.cbw,] 
   
   outcomes.impute.endog.cbw <- list("M"=best.var.outcome.cbw.m.imputed, # imputed data
                               "mask"=mask.cbw, 
                               "W"= p.weights.cbw.equal, # weights are equal 
                               "z.cbw.eastern"=z.cbw.eastern,
                               "z.cbw.swiss"= z.cbw.swiss,
                               "eastern"=eastern.cluster.cbw, 
                               "swiss"=swiss.cluster.cbw)
   
   impute.endog.best.var.outcome.cbw.eastern <- MCEst(outcomes.impute.endog.cbw, cluster='eastern', rev=TRUE, covars=FALSE)
   impute.endog.best.var.outcome.cbw.swiss <- MCEst(outcomes.impute.endog.cbw, cluster='swiss', rev=TRUE, covars=FALSE)
   best.var.outcome.cbw.hat <- matrix(NA, nrow(mask.cbw),ncol(mask.cbw),
                                            dimnames = list(rownames(mask.cbw), colnames(mask.cbw)))
   best.var.outcome.cbw.hat[!rownames(best.var.outcome.cbw.hat)%in% swiss.cluster.cbw,] <- best.var.outcome.cbw.m.imputed[!rownames(best.var.outcome.cbw.m.imputed)%in% swiss.cluster.cbw,]*(1-mask.cbw[!rownames(mask.cbw)%in% swiss.cluster.cbw,]) + 
     impute.endog.best.var.outcome.cbw.eastern$Mhat*mask.cbw[!rownames(mask.cbw)%in% swiss.cluster.cbw,]  # only endogenous values imputed
   
   best.var.outcome.cbw.hat[rownames(best.var.outcome.cbw.hat)%in% swiss.cluster.cbw,] <- (best.var.outcome.cbw.m.imputed[!rownames(best.var.outcome.cbw.m.imputed)%in% eastern.cluster.cbw,]*(1-mask.cbw[!rownames(mask.cbw)%in% eastern.cluster.cbw,]) + 
     impute.endog.best.var.outcome.cbw.swiss$Mhat*mask.cbw[!rownames(mask.cbw)%in% eastern.cluster.cbw,])[swiss.cluster.cbw,] 
   
   colnames(best.var.outcome.cbw.hat) <- colnames(mask.cbw)
   rownames(best.var.outcome.cbw.hat) <- rownames(mask.cbw)
   
  ## Estimate propensity scores by matrix completion
   
   propensity.model.cbw.data <- list("M"=mask.cbw, 
                                     "X" = best.var.outcome.cbw.m.imputed, # var with endogenous values + imputed missing values
                                     "X.hat"= best.var.outcome.cbw.hat, # var with imputed endogenous values
                                     "mask" = matrix(0, nrow(mask.cbw),ncol(mask.cbw),
                                                     dimnames = list(rownames(mask.cbw), colnames(mask.cbw))), # no missing entries
                                     "z.cbw.eastern"=z.cbw.eastern,
                                     "z.cbw.swiss"=z.cbw.swiss,
                                     "W"= p.weights.cbw.equal,
                                     "eastern"=eastern.cluster.cbw, 
                                     "swiss"=swiss.cluster.cbw) 
   
   propensity.model.cbw.eastern <- MCEst(propensity.model.cbw.data, cluster='eastern', rev=TRUE, covars=TRUE, prop.model=TRUE) 
   propensity.model.cbw.swiss <- MCEst(propensity.model.cbw.data, cluster='swiss', rev=TRUE, covars=TRUE, prop.model=TRUE) 
   
   propensity.model.cbw.values <- matrix(NA, nrow(mask.cbw),ncol(mask.cbw),
                                      dimnames = list(rownames(mask.cbw), colnames(mask.cbw)))
   
   propensity.model.cbw.values[!rownames(propensity.model.cbw.values)%in% swiss.cluster.cbw,] <- propensity.model.cbw.eastern$Mhat # L + X_hat*B + u + v
   propensity.model.cbw.values[rownames(propensity.model.cbw.values)%in% swiss.cluster.cbw,] <- propensity.model.cbw.swiss$Mhat[rownames(propensity.model.cbw.swiss$Mhat)%in% swiss.cluster.cbw,] # L + X_hat*B + u + v
   
   colnames(propensity.model.cbw.values) <- colnames(mask.cbw)
   rownames(propensity.model.cbw.values) <- rownames(mask.cbw)
  
   p.weights.cbw <- matrix(propensity.model.cbw.values, 
                          nrow = nrow(data.cbw ), 
                          ncol= ncol(data.cbw),
                          dimnames = list(rownames(data.cbw), colnames(data.cbw)), byrow = TRUE) # (N x T)
  # Save
  
  outcomes.cbw <- list("M"=data.cbw, "mask"=mask.cbw, 
                       "W"= p.weights.cbw, "W.equal"=p.weights.cbw.equal,
                       "z.cbw.eastern"= z.cbw.eastern, "z.cbw.swiss"= z.cbw.swiss,
                       "X"=best.var.outcome.cbw.m.imputed, "X.hat"=best.var.outcome.cbw.hat, "X.missing" = best.var.outcome.cbw.m,
                       "mc.outcome.eastern"=impute.best.var.outcome.cbw.eastern, "mc.outcome.swiss"=impute.best.var.outcome.cbw.swiss, 
                       "mc.propensity.eastern"=propensity.model.cbw.eastern, "mc.propensity.swiss"=propensity.model.cbw.swiss,
                               "treated"=rownames(mask.cbw)[rownames(mask.cbw)%in%switch.treated.cbw],
                               "control"=rownames(mask.cbw)[rownames(mask.cbw)%in%always.treated.cbw],
                               "eastern"=eastern.cluster.cbw, "swiss"=swiss.cluster.cbw)
  
  saveRDS(outcomes.cbw, paste0("data/outcomes-cbw-",o,".rds"))
}