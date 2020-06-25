## Load aggregated data

library(dplyr)
library(readstata13)
library(glmnet)
library(caret)

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

covariates <- c("AV_age_22","AV_age_27","AV_age_32","AV_age_37","AV_age_42","AV_age_47","AV_age_52","AV_age_57","AV_women","AV_lowEDU","AV_mediumEDU","AV_highEDU","AV_migr","HHincome_COUNTRY","Lang",               
                "XCGvsEURO_05","GDPcapitaR","R_GDPcapitaR","density","R_HHincome_COUNTRY2","R_HHincome_COUNTRY","AV_HHsize","AV_HHsizeCHILD","AV_HHsizeOLD","AV_single","AV_oneadultNOCHILD","AV_oneadultCHILD",   
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

covars <- lapply(covariates, function(c){
  mat <- reshape(data.frame(data[c("REGION","yearquarter",c)]), idvar = "REGION", timevar = "yearquarter", direction = "wide") # N x T
  colnames(mat) <- sub(paste0(c,"."),"", colnames(mat))
  rownames(mat) <- mat$REGION
  mat <- mat[,-1]
})

# Time-invariant covariates

covars.cbw <- as.matrix(bind_cols(sapply(1:length(covars), function(i){ # only post-treatment years
  covars[[i]][,which(colnames(covars[[i]])=='20111'):ncol(covars[[i]])][rownames(covars[[i]]) %in% cbw,]
}))) # N x # predictors
rownames(covars.cbw) <- rownames(covars[[1]][rownames(covars[[1]]) %in% cbw,])

covars.lm <- as.matrix(bind_cols(sapply(1:length(covars), function(i){ # only pre-treatment years
  covars[[i]][,1:which(colnames(covars[[i]])=='20071')][rownames(covars[[i]]) %in% lm,]
}))) # N x # predictors
rownames(covars.lm) <- rownames(covars[[1]][rownames(covars[[1]]) %in% lm,])

# Remove duplicated columns

covars.cbw <- covars.cbw[,!duplicated(t(covars.cbw))]

covars.lm <- covars.lm[,!duplicated(t(covars.lm))]

# Impute missing with column medians

covars.cbw <- predict(preProcess(covars.cbw, method = c("medianImpute")), covars.cbw)

covars.lm <- predict(preProcess(covars.lm, method = c("medianImpute")), covars.lm)

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
  
  ## estimate propensity scores with time-invariant covariates
  
  # CBW
  
  covars.cbw.reduced <- covars.cbw[rownames(covars.cbw)%in%rownames(mask.cbw),]

  logitMod.cbw <- cv.glmnet(x=covars.cbw.reduced,
                            y=as.factor(mask.cbw[,"20084"]), # 1 if switch.treated.cbw
                            family="binomial", nfolds= 10, parallel = TRUE)

  covars.cbw.preds <- as.vector(predict(logitMod.cbw, covars.cbw.reduced, type="response", s ="lambda.min"))

  names(covars.cbw.preds) <- names(as.factor(mask.cbw[,"20084"]))
  
  z.cbw.eastern <- round(c(seq(1, 0.7, length.out=which(colnames(mask.cbw)=="20111")),
                     seq(0.712, 1, length.out=ncol(mask.cbw)-which(colnames(mask.cbw)=="20111"))),3) # elapsed time since treatment
  
  z.cbw.swiss <- round(c(seq(1, 0.7, length.out=which(colnames(mask.cbw)=="20091")),
                           seq(0.719, 1, length.out=ncol(mask.cbw)-which(colnames(mask.cbw)=="20091"))),3)
  
  p.weights.cbw <- matrix(0, nrow = nrow(data.cbw ), 
                                        ncol= ncol(data.cbw),
                                        dimnames = list(rownames(data.cbw), colnames(data.cbw))) # (N x T)
  
  p.weights.cbw[rownames(p.weights.cbw) %in% eastern.cluster.cbw,] <- covars.cbw.preds[names(covars.cbw.preds) %in% eastern.cluster.cbw]%*%t(z.cbw.eastern) # inner product
  p.weights.cbw[rownames(p.weights.cbw) %in% swiss.cluster.cbw,] <- covars.cbw.preds[names(covars.cbw.preds) %in% swiss.cluster.cbw]%*%t(z.cbw.swiss) # inner product
  p.weights.cbw[rownames(p.weights.cbw) %in% always.treated.cbw,] <- covars.cbw.preds[names(covars.cbw.preds) %in% always.treated.cbw] # no time adjustment for controls
  
  # LM
  
  covars.lm.reduced <- covars.lm[rownames(covars.lm)%in%rownames(mask.lm),]
  
  logitMod.lm <- cv.glmnet(x=covars.lm.reduced,
                            y=as.factor(mask.lm[,"20111"]), # 1 if switch.treated.lm
                            family="binomial", nfolds= 10, parallel = TRUE)
  
  covars.lm.preds <- as.vector(predict(logitMod.lm, covars.lm.reduced, type="response", s ="lambda.min"))
  
  names(covars.lm.preds) <- names(as.factor(mask.lm[,"20111"]))
  
  z.lm.eastern <- round(c(seq(1, 0.7, length.out=which(colnames(mask.lm)=="20111")),
                           seq(0.712, 1, length.out=ncol(mask.lm)-which(colnames(mask.lm)=="20111"))),3) # elapsed time since treatment
  
  z.lm.swiss <- round(c(seq(1, 0.7, length.out=which(colnames(mask.lm)=="20091")),
                         seq(0.719, 1, length.out=ncol(mask.lm)-which(colnames(mask.lm)=="20091"))),3)
  
  p.weights.lm <- matrix(0, nrow = nrow(data.lm ), 
                          ncol= ncol(data.lm),
                          dimnames = list(rownames(data.lm), colnames(data.lm))) # (N x T)
  
  p.weights.lm[rownames(p.weights.lm) %in% eastern.cluster.lm,] <- covars.lm.preds[names(covars.lm.preds) %in% eastern.cluster.lm]%*%t(z.lm.eastern) # inner product
  p.weights.lm[rownames(p.weights.lm) %in% swiss.cluster.lm,] <- covars.lm.preds[names(covars.lm.preds) %in% swiss.cluster.lm]%*%t(z.lm.swiss) # inner product
  p.weights.lm[rownames(p.weights.lm) %in% never.treated.lm,] <- covars.lm.preds[names(covars.lm.preds) %in% never.treated.lm] # no time adjustment for controls
  
  # Save
  
  outcomes.cbw <- list("M"=data.cbw, "mask"=mask.cbw, "W"= p.weights.cbw, "X"=covars.cbw.reduced, 
                               "treated"=rownames(mask.cbw)[rownames(mask.cbw)%in%switch.treated.cbw],
                               "control"=rownames(mask.cbw)[rownames(mask.cbw)%in%always.treated.cbw])
  outcomes.lm <- list("M"=data.lm, "mask"=mask.lm, "W"= p.weights.lm, "X"=covars.lm.reduced, 
                              "treated"=rownames(mask.lm)[rownames(mask.lm)%in%switch.treated.lm],
                              "control"=rownames(mask.lm)[rownames(mask.lm)%in%never.treated.lm])
  
  saveRDS(outcomes.cbw, paste0("data/outcomes-cbw-",o,".rds"))
  saveRDS(outcomes.lm, paste0("data/outcomes-lm-",o,".rds"))
}