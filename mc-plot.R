######################################################################
# Plot Time-series and causal impacts, and estimated latent trends #
######################################################################

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
source("getTrends.R")
source("TsPlotTrends.R")

PlotMCCapacity <- function(observed,main,y.title,mc_est_eastern,mc_est_swiss,boot_result_eastern,boot_result_swiss,treated,control,eastern,swiss,vline,vline2,breaks,labels,att.label,censor=FALSE){
  ## Create Timeseries data
  
  predicted.eastern <- mc_est_eastern$Mhat
  pointwise.eastern <- mc_est_eastern$impact # = boot_result$t0
  
  predicted.swiss <- mc_est_swiss$Mhat
  pointwise.swiss <- mc_est_swiss$impact # = boot_result$t0
  
  # Create confidence intervals
  
  boot.ci.lower.eastern <-matrix(apply(boot_result_eastern$t, 2, function(x) quantile(x, 0.025)), nrow=dim(pointwise.eastern)[1], ncol=dim(pointwise.eastern)[2], byrow=FALSE,
                              dimnames=dimnames(pointwise.eastern)) # alpha/2
  boot.ci.upper.eastern <-matrix(apply(boot_result_eastern$t, 2, function(x) quantile(x, 0.975)), nrow=dim(pointwise.eastern)[1], ncol=dim(pointwise.eastern)[2], byrow=FALSE,
                              dimnames=dimnames(pointwise.eastern)) # 1-alpha/2
  
  boot.ci.lower.swiss <-matrix(apply(boot_result_swiss$t, 2, function(x) quantile(x, 0.025)), nrow=dim(pointwise.swiss)[1], ncol=dim(pointwise.swiss)[2], byrow=FALSE,
                                 dimnames=dimnames(pointwise.swiss)) # alpha/2
  boot.ci.upper.swiss <-matrix(apply(boot_result_swiss$t, 2, function(x) quantile(x, 0.975)), nrow=dim(pointwise.swiss)[1], ncol=dim(pointwise.swiss)[2], byrow=FALSE,
                                 dimnames=dimnames(pointwise.swiss)) # 1-alpha/2
  
  ## Plot Timeseries 
  
  treat.status <- matrix(rownames(observed), nrow=nrow(observed), ncol=1)
  treat.status[rownames(observed) %in% eastern] <- "eastern"
  treat.status[rownames(observed) %in% swiss] <- "swiss"
  treat.status[rownames(observed) %in% control] <- "control"
  treat.status <- matrix(treat.status, dimnames=list(NULL, "status"))
  
  treat.status.eastern <- matrix(rownames(pointwise.eastern), nrow=nrow(pointwise.eastern), ncol=1)
  treat.status.eastern[rownames(pointwise.eastern) %in% eastern] <- "eastern"
  treat.status.eastern[rownames(pointwise.eastern) %in% control] <- "control"
  treat.status.eastern <- matrix(treat.status.eastern, dimnames=list(NULL, "status"))
  
  treat.status.swiss <- matrix(rownames(pointwise.swiss), nrow=nrow(pointwise.swiss), ncol=1)
  treat.status.swiss[rownames(pointwise.swiss) %in% swiss] <- "swiss"
  treat.status.swiss[rownames(pointwise.swiss) %in% control] <- "control"
  treat.status.swiss <- matrix(treat.status.swiss, dimnames=list(NULL, "status"))
  
  observed.mean <-  aggregate(observed, list(treat.status), mean)[-1]
  
  predicted.mean.eastern <-  aggregate(predicted.eastern, list(treat.status.eastern), mean)[-1]
  pointwise.mean.eastern <- aggregate(pointwise.eastern, list(treat.status.eastern), mean)[-1]
  pointwise.ci.lower.mean.eastern <- aggregate(boot.ci.lower.eastern, list(treat.status.eastern), mean)[-1]
  pointwise.ci.upper.mean.eastern <- aggregate(boot.ci.upper.eastern, list(treat.status.eastern), mean)[-1]
  
  predicted.mean.swiss <-  aggregate(predicted.swiss, list(treat.status.swiss), mean)[2,][-1]
  pointwise.mean.swiss <- aggregate(pointwise.swiss, list(treat.status.swiss), mean)[2,][-1]
  pointwise.ci.lower.mean.swiss <- aggregate(boot.ci.lower.swiss, list(treat.status.swiss), mean)[2,][-1]
  pointwise.ci.upper.mean.swiss <- aggregate(boot.ci.upper.swiss, list(treat.status.swiss), mean)[2,][-1]

  ts.means <- cbind(t(observed.mean), t(predicted.mean.eastern), t(predicted.mean.swiss), t(pointwise.mean.eastern), t(pointwise.mean.swiss))
  colnames(ts.means) <- c("observed.control","observed.eastern","observed.swiss","predicted.control","predicted.eastern","predicted.swiss","pointwise.control","pointwise.eastern","pointwise.swiss")
  ts.means <- cbind(ts.means, "year"=as.numeric(rownames(ts.means)))
  ts.means.m <- melt(data.frame(ts.means), id.var=c("year"))
  
  ts.ci.lower.means <- cbind(t(pointwise.ci.lower.mean.eastern),t(pointwise.ci.lower.mean.swiss))
  colnames(ts.ci.lower.means) <- c("pointwise.control","pointwise.eastern","pointwise.swiss")
  ts.ci.lower.means <- cbind(ts.ci.lower.means, "year"=as.numeric(rownames(ts.means)))
  ts.ci.lower.means.m <- melt(data.frame(ts.ci.lower.means), id.var=c("year"))
  
  ts.ci.upper.means <- cbind(t(pointwise.ci.upper.mean.eastern),t(pointwise.ci.upper.mean.swiss))
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
  
  if(censor){
    ts.means.m$value[ts.means.m$year > vline2 & (ts.means.m$variable=="pointwise.eastern" | ts.means.m$variable=="predicted.eastern")] <- NA # censor
    ts.means.m$value[ts.means.m$year > vline & (ts.means.m$variable=="pointwise.swiss" | ts.means.m$variable=="predicted.swiss")] <- NA

    ts.means.m$upper[ts.means.m$year > vline2 & (ts.means.m$variable=="pointwise.eastern" | ts.means.m$variable=="predicted.eastern")] <- NA
    ts.means.m$upper[ts.means.m$year > vline & (ts.means.m$variable=="pointwise.swiss" | ts.means.m$variable=="predicted.swiss")] <- NA

    ts.means.m$lower[ts.means.m$year > vline2 & ts.means.m$variable=="pointwise.eastern"] <- NA
    ts.means.m$lower[ts.means.m$year > vline & ts.means.m$variable=="pointwise.swiss"] <- NA
  }
  
  ts.means.m$quarter <- rep(1:nrow(ts.means), each=9) # for x axis
  
  ts.plot <- TsPlot(df=ts.means.m,main=main, y.title=y.title,vline,vline2,breaks,labels,hline=ts.means.m$hline)
  
  return(ts.plot)
}

## Plot time-series

outcome.vars <- c("CBWbord","CBWbordEMPL")
outcomes.labels <- c("% working in border region",
                     "% working in border region, conditional on employment")

covarflag <- c("-covars")

for(o in outcome.vars){
  for(cf in covarflag){
    
    outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
    
    mc.estimates.cbw.eastern <- readRDS(paste0("results/mc-estimates-cbw-eastern-",o,cf,".rds"))
    boot.cbw.eastern <- readRDS(paste0("results/boot-cbw-eastern-",o,cf,".rds"))
    
    mc.estimates.cbw.swiss <- readRDS(paste0("results/mc-estimates-cbw-swiss-",o,cf,".rds"))
    boot.cbw.swiss <- readRDS(paste0("results/boot-cbw-swiss-",o,cf,".rds"))
    
    mc.plot <- PlotMCCapacity(observed = outcomes.cbw$M, 
                              y.title=outcomes.labels[which(outcome.vars==o)],
                           #   main = "Retrospective prediction for later-treated, by cluster",
                              main= "",
                              mc_est_eastern=mc.estimates.cbw.eastern, 
                              mc_est_swiss=mc.estimates.cbw.swiss,
                              boot_result_eastern=boot.cbw.eastern,
                              boot_result_swiss=boot.cbw.swiss,
                              treated=outcomes.cbw$treated, 
                              control=outcomes.cbw$control, 
                              eastern=outcomes.cbw$eastern,
                              swiss= outcomes.cbw$swiss,
                              vline=16,vline2=24,
                              breaks=c(4,12,20,28,36,44,52,60),
                              labels=c("2005Q4","2007Q4","2009Q4","2011Q4","2013Q4","2015Q4","2017Q4","2019Q4"),
                              att.label = TeX("$\\hat{\\tau}_{t}^{LT}$"))
    
    ggsave(paste0("plots/mc-estimates-cbw-",o,cf,".png"), mc.plot + theme(legend.position = "none"), scale=1.25)
    ggsave(paste0("plots/mc-estimates-cbw-",o,cf,"-slides.png"), mc.plot + ggtitle("Matrix completion estimates of the effect of Schengen and FoM") + theme(legend.position = "none", plot.title = element_text(family="serif", size=16, hjust = 0.5)), scale=1.25) 
  
    # Plot estimated latent trends of the Swiss and Eastern blocks
    
    trend.eastern <- getTrends(L=mc.estimates.cbw.eastern$L)
    trend.swiss <- getTrends(L=mc.estimates.cbw.swiss$L)
    
    #Plot first trend
    trend.first.data <- data.frame("eastern"=trend.eastern$first,
                             "swiss"=trend.swiss$first,
                               "year"=colnames(outcomes.cbw$M))
    
    trend.first.data.m <- melt(trend.first.data,id="year")
    
    trend.first.data.m$quarter <- rep(1:nrow(trend.first.data), times=2) # for x axis
      
    trend.first.plot <- TsPlotTrends(df=trend.first.data.m,main="", 
                               y.title=outcomes.labels[which(outcome.vars==o)],
                               vline=16,
                               vline2=24,
                               breaks=c(4,12,20,28,36,44,52,60),
                               labels=c("2005Q4","2007Q4","2009Q4","2011Q4","2013Q4","2015Q4","2017Q4","2019Q4"))
    
    
    ggsave(paste0("plots/mc-trend-first-",o,cf,".png"), trend.first.plot + theme(legend.position = "none"), scale=1.25)
    ggsave(paste0("plots/mc-trend-first-",o,cf,"-slides.png"), trend.first.plot + ggtitle("Matrix completion estimates: first latend trend") + theme(legend.position = "none", plot.title = element_text(family="serif", size=16, hjust = 0.5)), scale=1.25) 
    
    #Plot second trend
    trend.second.data <- data.frame("eastern"=trend.eastern$second,
                                   "swiss"=trend.swiss$second,
                                   "year"=colnames(outcomes.cbw$M))
    
    trend.second.data.m <- melt(trend.second.data,id="year")
    
    trend.second.data.m$quarter <- rep(1:nrow(trend.second.data), times=2) # for x axis
    
    trend.second.plot <- TsPlotTrends(df=trend.second.data.m,main="", 
                                     y.title=outcomes.labels[which(outcome.vars==o)],
                                     vline=16,
                                     vline2=24,
                                     breaks=c(4,12,20,28,36,44,52,60),
                                     labels=c("2005Q4","2007Q4","2009Q4","2011Q4","2013Q4","2015Q4","2017Q4","2019Q4"))
    
    
    ggsave(paste0("plots/mc-trend-second-",o,cf,".png"), trend.second.plot + theme(legend.position = "none"), scale=1.25)
    ggsave(paste0("plots/mc-trend-second-",o,cf,"-slides.png"), trend.second.plot + ggtitle("Matrix completion estimates: second latent trend") + theme(legend.position = "none", plot.title = element_text(family="serif", size=16, hjust = 0.5)), scale=1.25) 
    }
} 