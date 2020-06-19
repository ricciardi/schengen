###################################
# Plot time-series and causal impacts #
###################################

require(reshape2)
require(dplyr)
require(zoo)
require(matrixStats)
require(tseries)
require(ggplot2)
library(latex2exp)

source("TsPlot.R")

PlotMCCapacity <- function(observed,main,y.title,t0,mc_est,boot_result,treated,control,vline,breaks,labels,att.label,rev=TRUE){
  ## Create time series data
  
  predicted <- mc_est$Mhat
  
  if(rev){
    pointwise <- abs(mc_est$impact)
  }
  else{
    pointwise <- mc_est$impact # boot_result$t0
  }
  
  pointwise.se <- matrix(apply(boot_result$t, 2, sd), nrow=dim(pointwise)[1], ncol=dim(pointwise)[2], byrow=FALSE)

  m <- ncol(observed)
  n <- t0
  
  ## Plot time series 
  
  treat.status <- matrix(rownames(pointwise), nrow=nrow(pointwise), ncol=1)
  treat.status[rownames(pointwise) %in% treated] <- "treated"
  treat.status[rownames(pointwise) %in% control] <- "control"
  treat.status <- matrix(treat.status, dimnames=list(NULL, "status"))
  
  observed.mean <-  aggregate(observed, list(treat.status), mean)[-1]
  predicted.mean <-  aggregate(predicted, list(treat.status), mean)[-1]
  pointwise.mean <- aggregate(pointwise, list(treat.status), mean, na.rm=TRUE)[-1]
  pointwise.se.mean <- aggregate(pointwise.se, list(treat.status), mean)[-1]
  
  ts.means <- cbind(t(observed.mean), t(predicted.mean), t(pointwise.mean))
  colnames(ts.means) <- c("observed.control","observed.treated","predicted.control","predicted.treated","pointwise.control","pointwise.treated")
  ts.means <- cbind(ts.means, "year"=as.numeric(rownames(ts.means)))
  ts.means.m <- melt(data.frame(ts.means), id.var=c("year"))
  
  ts.se.means <- cbind(t(pointwise.se.mean))
  colnames(ts.se.means) <- c("pointwise.control","pointwise.treated")
  ts.se.means <- cbind(ts.se.means, "year"=as.numeric(rownames(ts.means)))
  ts.se.means.m <- melt(data.frame(ts.se.means), id.var=c("year"))
  
  ts.means.m <- merge(ts.means.m, ts.se.means.m, by=c("year","variable"), all.x=TRUE) # bind std. error
  colnames(ts.means.m) <- c("year", "variable", "value", "se")
  
  ts.means.m <- ts.means.m %>%
    mutate(upper = value + 1.96*se,
           lower = value - 1.96*se)
  
  # Labels
  
  ts.means.m$series <- NA
  ts.means.m$series[grep("observed.", ts.means.m$variable)] <- "Time-series"
  ts.means.m$series[grep("predicted.", ts.means.m$variable)] <- "Time-series"
  ts.means.m$series[grep("pointwise.", ts.means.m$variable)] <- att.label
  
  ts.means.m$series<- factor(ts.means.m$series, levels=c("Time-series", att.label)) # reverse order
  
  ts.plot <- TsPlot(df=ts.means.m,main=main, y.title=y.title,vline,breaks,labels)
  
  return(ts.plot)
}

outcomes <- c("CBWbord","CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_0","seekdur_1_2","seekdur_3more")
outcomes.labels <- c()

for(o in outcomes){
  print(o)
  
  ## Analysis 1: ST vs AT (retrospective, X=CBW) 
  
  # Eastern cluster
  
  outcomes.cbw.eastern <- readRDS(paste0("data/outcomes-cbw-eastern-",o,".rds"))
  mc.estimates.cbw.eastern <- readRDS(paste0("results/mc-estimates-cbw-eastern-",o,".rds"))
  boot.cbw.eastern <- readRDS(paste0("results/boot-cbw-eastern-",o,".rds"))
  
  mc.plot <- PlotMCCapacity(observed = outcomes.east$M, 
                            y.title=outcomes.labels[which(outcomes==o)],
                            main = "Retrospective prediction for switch-treated",
                            t0=which(colnames(outcomes.east$M)=="20081"),
                            mc_est=MC_estimates, boot_result=boot, treated=outcomes.east$treated, control=outcomes.east$control, vline=20081,
                            breaks=c(20042,20081,20121,20161,20184),
                            labels=c(20042,20081,20121,20161,20184),
                            att.label = TeX("$\\tau = E(\\hat{Y}(1)_{LT, pre} - Y(0)_{LT,pre})$"),
                            rev=TRUE)
  
  ggsave("plots/mc.png", mc.plot, width=8.5, height=11)
}