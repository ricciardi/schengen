## Load aggregated data

library(dplyr)
library(readstata13)

data <- read.dta13("../FINAL.dta", generate.factors=T) # includes variables for 2 analyses + covariates 

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
switch.treated.lm <- unique(data$REGION[which(data$treated_LM==c("2008: Schengen; 2011: FoM", "2007: FoM; 2009: Schengen"))])
length(never.treated.lm[!never.treated.lm%in%switch.treated.lm]) == length(never.treated.lm) # ensure these groups don't overlap

# Clusters
## Eastern cluster (treated_X=1)
### SCHENGEN_X=1 from 20081 and FoM_X=1 from 20111
## Swiss cluster (treated_X=2)
### FoM_X=1 from 2007Q2 and SCHENGEN_X=1 from 2009Q1

# no. results = 2*2*2*9 = 72

# create yearquarter
data$quarter <-NA
data$quarter <- rep(1:4, nrow(data)/4)
data$yearquarter <- as.numeric(paste0(data$year,data$quarter))

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
  mask.cbw.eastern <- matrix(0, nrow = nrow(data.cbw ), 
                             ncol= ncol(data.cbw),
                             dimnames = list(rownames(data.cbw), colnames(data.cbw))) # (N x T)
  
  mask.cbw.swiss <- matrix(0, nrow = nrow(data.cbw ), 
                           ncol= ncol(data.cbw),
                           dimnames = list(rownames(data.cbw), colnames(data.cbw))) # (N x T)
  
  mask.lm.eastern <- matrix(0, nrow = nrow(data.lm ), 
                            ncol= ncol(data.lm),
                            dimnames = list(rownames(data.lm), colnames(data.lm))) # (N x T)
  
  mask.lm.swiss <- matrix(0, nrow = nrow(data.lm ), 
                          ncol= ncol(data.lm),
                          dimnames = list(rownames(data.lm), colnames(data.lm))) # (N x T)
  
  for(i in switch.treated.cbw){
    mask.cbw.eastern[,colnames(mask.cbw.eastern)%in%data$yearquarter[data$REGION==i & data$SCHENGEN_CBW==1 & data$treated_CBW=="2008: Schengen; 2011: FoM"]][rownames(mask.cbw.eastern)%in%c(i),] <- 1 # SCHENGEN_X=1 from 20081 and FoM_X=1 from 20111
  }
  
  for(i in switch.treated.cbw){
    mask.cbw.swiss[,colnames(mask.cbw.swiss)%in%data$yearquarter[data$REGION==i & data$FoM_CBW==1 & data$treated_CBW=="2007: FoM; 2009: Schengen"]][rownames(mask.cbw.swiss)%in%c(i),] <- 1 # FoM_X=1 from 2007Q2 and SCHENGEN_X=1 from 2009Q1
  }
  
  for(i in switch.treated.lm){
    mask.lm.eastern[,colnames(mask.lm.eastern)%in%data$yearquarter[data$REGION==i & data$SCHENGEN_LM==1 & data$treated_LM=="2008: Schengen; 2011: FoM"]][rownames(mask.lm.eastern)%in%c(i),] <- 1 # SCHENGEN_X=1 from 20081 and FoM_X=1 from 20111
  }
  
  for(i in switch.treated.lm){
    mask.lm.swiss[,colnames(mask.lm.swiss)%in%data$yearquarter[data$REGION==i & data$FoM_LM==1 & data$treated_LM=="2007: FoM; 2009: Schengen"]][rownames(mask.lm.swiss)%in%c(i),] <- 1 # FoM_X=1 from 2007Q2 and SCHENGEN_X=1 from 2009Q1
  }
  # Save
  
  outcomes.cbw.eastern <- list("M"=data.cbw, "mask"=mask.cbw.eastern,"treated"=switch.treated.cbw,"control"=always.treated.cbw)
  outcomes.lm.eastern <- list("M"=data.lm, "mask"=mask.lm.eastern,"treated"=switch.treated.lm,"control"=never.treated.lm)
  
  outcomes.cbw.swiss <- list("M"=data.cbw, "mask"=mask.cbw.swiss,"treated"=switch.treated.cbw,"control"=always.treated.cbw)
  outcomes.lm.swiss <- list("M"=data.lm, "mask"=mask.lm.swiss,"treated"=switch.treated.lm,"control"=never.treated.lm)
  
  saveRDS(outcomes.cbw.eastern, paste0("data/outcomes-cbw-eastern-",o,".rds"))
  saveRDS(outcomes.lm.eastern, paste0("data/outcomes-lm-eastern-",o,".rds"))
  
  saveRDS(outcomes.cbw.swiss, paste0("data/outcomes-cbw-swiss-",o,".rds"))
  saveRDS(outcomes.lm.swiss,paste0("data/outcomes-lm-swiss-",o,".rds"))
}

# Covariates matrices

covars <- lapply(covariates, function(c){
  mat <- reshape(data.frame(data[c("REGION","yearquarter",c)]), idvar = "REGION", timevar = "yearquarter", direction = "wide") # N x T
  colnames(mat) <- sub(paste0(c,"."),"", colnames(mat))
  rownames(mat) <- mat$REGION
  mat <- mat[,-1]
})

names(covars) <- covariates
saveRDS(covars,"data/covariates.rds")