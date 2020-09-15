###################################
# Plot Time-series and causal impacts #
###################################

require(reshape2)
require(dplyr)
require(zoo)
require(matrixStats)
require(tseries)
require(ggplot2)
library(latex2exp)
library(wesanderson)
library(boot)

source("TsPlot.R")

PlotMCCapacity <- function(observed,main,y.title,mc_est,boot_result,treated,control,eastern,swiss,vline,vline2,breaks,labels,att.label,rev){
  ## Create time series data
  
  predicted <- mc_est$Mhat
  
  pointwise <- mc_est$impact # = boot_result$t0
  
  # Create confidence intervals
  
  boot.ci.lower <-matrix(apply(boot_result$t, 2, function(x) quantile(x, 0.025)), nrow=dim(pointwise)[1], ncol=dim(pointwise)[2], byrow=FALSE,
                              dimnames=dimnames(pointwise)) # alpha/2
  boot.ci.upper <-matrix(apply(boot_result$t, 2, function(x) quantile(x, 0.975)), nrow=dim(pointwise)[1], ncol=dim(pointwise)[2], byrow=FALSE,
                              dimnames=dimnames(pointwise)) # 1-alpha/2
  
  ## Plot time series 
  
  treat.status <- matrix(rownames(pointwise), nrow=nrow(pointwise), ncol=1)
  treat.status[rownames(pointwise) %in% eastern] <- "eastern"
  treat.status[rownames(pointwise) %in% swiss] <- "swiss"
  treat.status[rownames(pointwise) %in% control] <- "control"
  treat.status <- matrix(treat.status, dimnames=list(NULL, "status"))
  
  observed.mean <-  aggregate(observed, list(treat.status), mean)[-1]
  predicted.mean <-  aggregate(predicted, list(treat.status), mean)[-1]
  pointwise.mean <- aggregate(pointwise, list(treat.status), mean)[-1]
  pointwise.ci.lower.mean <- aggregate(boot.ci.lower, list(treat.status), mean)[-1]
  pointwise.ci.upper.mean <- aggregate(boot.ci.upper, list(treat.status), mean)[-1]

  ts.means <- cbind(t(observed.mean), t(predicted.mean), t(pointwise.mean))
  colnames(ts.means) <- c("observed.control","observed.eastern","observed.swiss","predicted.control","predicted.eastern","predicted.swiss","pointwise.control","pointwise.eastern","pointwise.swiss")
  ts.means <- cbind(ts.means, "year"=as.numeric(rownames(ts.means)))
  ts.means.m <- melt(data.frame(ts.means), id.var=c("year"))
  
  ts.ci.lower.means <- t(pointwise.ci.lower.mean)
  colnames(ts.ci.lower.means) <- c("pointwise.control","pointwise.eastern","pointwise.swiss")
  ts.ci.lower.means <- cbind(ts.ci.lower.means, "year"=as.numeric(rownames(ts.means)))
  ts.ci.lower.means.m <- melt(data.frame(ts.ci.lower.means), id.var=c("year"))
  
  ts.ci.upper.means <- t(pointwise.ci.upper.mean)
  colnames(ts.ci.upper.means) <- c("pointwise.control","pointwise.eastern","pointwise.swiss")
  ts.ci.upper.means <- cbind(ts.ci.upper.means, "year"=as.numeric(rownames(ts.means)))
  ts.ci.upper.means.m <- melt(data.frame(ts.ci.upper.means), id.var=c("year"))
  
  ts.means.m <- merge(ts.means.m, ts.ci.lower.means.m, by=c("year","variable"), all.x=TRUE) # bind std. error
  colnames(ts.means.m) <- c("year", "variable", "value", "boot.lower")

  ts.means.m <- merge(ts.means.m, ts.ci.upper.means.m, by=c("year","variable"), all.x=TRUE) # bind std. error
  colnames(ts.means.m) <- c("year", "variable", "value", "boot.lower","boot.upper")
  
  ts.means.m <- ts.means.m %>%
    mutate(upper = value -boot.lower,
           lower = value -boot.upper)
  
  # Labels
  
  ts.means.m$series <- NA
  ts.means.m$series[grep("observed.", ts.means.m$variable)] <- "Timeseries"
  ts.means.m$series[grep("predicted.", ts.means.m$variable)] <- "Timeseries"
  ts.means.m$series[grep("pointwise.", ts.means.m$variable)] <- att.label
  
  ts.means.m$series<- factor(ts.means.m$series, levels=c("Timeseries", att.label)) # reverse order
  
  ts.means.m$hline <-NA
  ts.means.m$hline[ts.means.m$series!="Timeseries"] <-0
  
  if(rev){
    ts.means.m$value[ts.means.m$year > vline2 & (ts.means.m$variable=="pointwise.eastern" | ts.means.m$variable=="predicted.eastern")] <- NA # censor 
    ts.means.m$value[ts.means.m$year > vline & (ts.means.m$variable=="pointwise.swiss" | ts.means.m$variable=="predicted.swiss")] <- NA
    
    ts.means.m$upper[ts.means.m$year > vline2 & (ts.means.m$variable=="pointwise.eastern" | ts.means.m$variable=="predicted.eastern")] <- NA 
    ts.means.m$upper[ts.means.m$year > vline & (ts.means.m$variable=="pointwise.swiss" | ts.means.m$variable=="predicted.swiss")] <- NA
    
    ts.means.m$lower[ts.means.m$year > vline2 & ts.means.m$variable=="pointwise.eastern"] <- NA 
    ts.means.m$lower[ts.means.m$year > vline & ts.means.m$variable=="pointwise.swiss"] <- NA
  }
  
  ts.plot <- TsPlot(df=ts.means.m,main=main, y.title=y.title,vline,vline2,breaks,labels,hline=ts.means.m$hline,rev)
  
  return(ts.plot)
}

## Plot time-series

outcome.vars <- c("CBWbord","CBWbordEMPL","Thwusual")
outcomes.labels <- c("% working in border region",
                     "% working in border region,\n conditional on employment",
                     "Average total working hours")

covarflag <- c("","-covars")

for(o in outcome.vars){
  for(cf in covarflag){
    
    ## Analysis 1: ST vs AT (retrospective, X=CBW) 
    
    outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
    mc.estimates.cbw <- readRDS(paste0("results/mc-estimates-cbw-",o,cf,".rds"))
    
    boot.cbw <- readRDS(paste0("results/boot-cbw-",o,cf,".rds"))
    
    mc.plot <- PlotMCCapacity(observed = outcomes.cbw$M, 
                              y.title=outcomes.labels[which(outcome.vars==o)],
                           #   main = "Retrospective prediction for later-treated, by cluster",
                              main= "",
                              mc_est=mc.estimates.cbw, 
                              boot_result=boot.cbw,
                              treated=outcomes.cbw$treated, 
                              control=outcomes.cbw$control, 
                              eastern =outcomes.cbw$eastern,
                              swiss= outcomes.cbw$swiss,
                              vline=20091,vline2=20111,
                              breaks=c(20051,20072,20091,20111,20184),
                              labels=c("2005Q1","2007Q2","2009Q1","2011Q1","2018Q4"),
                              att.label = TeX("$\\hat{\\tau}_{t}^{\ATT}$"),
                              rev=TRUE)
    
    ggsave(paste0("plots/mc-estimates-cbw-",o,cf,".png"), mc.plot, width=8.5, height=11)
    ggsave(paste0("plots/mc-estimates-cbw-",o,cf,"-slides.png"), mc.plot + ggtitle("Matrix completion estimates of the effect of Schengen and FoM") + theme(plot.title = element_text(family="serif", size=16, hjust = 0.5)), scale=2) 
  }
} 