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

outcome.vars <- c("CBWbordEMPL","empl","Thwusual","unempl","seekdur_3more")
outcomes.labels <- c("% working in border region",
                     "Employment rate",
                     "Average total working hours",
                     "Unemployment rate",
                     "% unemployed for < 1 year")

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
                              att.label = TeX("$\\hat{\\bar{\\tau}}_t$"),
                              rev=TRUE)
    
    ggsave(paste0("plots/mc-estimates-cbw-",o,cf,".png"), mc.plot, width=8.5, height=11)
    ggsave(paste0("plots/mc-estimates-cbw-",o,cf,"-slides.png"), mc.plot, scale=1.75) 
  }
} 

# ## Plot trajectory CIs
# 
# ## Analysis 1: ST vs AT (retrospective, X=CBW) 
# 
# boot.trajectory.eastern.placebo.cbw  <- lapply(outcome.vars, function(o,cf="-covars"){
#   boot  <- readRDS(paste0("results/boot-cbw-trajectory-eastern-",o,cf,".rds"))
#   return(list("t0"=boot$t0, "ci.lower"=boot.ci(boot, type="perc")$percent[4],
#               "ci.upper"=boot.ci(boot, type="perc")$percent[5]))
# })
# names(boot.trajectory.eastern.placebo.cbw ) <- outcome.vars
# 
# boot.trajectory.swiss.placebo.cbw  <- lapply(outcome.vars, function(o,cf="-covars"){
#   boot  <- readRDS(paste0("results/boot-cbw-trajectory-swiss-",o,cf,".rds"))
#   return(list("t0"=boot$t0, "ci.lower"=boot.ci(boot, type="perc")$percent[4],
#               "ci.upper"=boot.ci(boot, type="perc")$percent[5]))
# })
# names(boot.trajectory.swiss.placebo.cbw ) <- outcome.vars
# 
# ci.values <- data.frame("t0"=c(sapply(outcome.vars, function(i) boot.trajectory.eastern.placebo.cbw[[i]]$t0),
#                                sapply(outcome.vars, function(i) boot.trajectory.swiss.placebo.cbw[[i]]$t0)),
#                         "lower"=c(sapply(outcome.vars, function(i) boot.trajectory.eastern.placebo.cbw[[i]]$ci.lower),
#                                   sapply(outcome.vars, function(i) boot.trajectory.swiss.placebo.cbw[[i]]$ci.lower)),
#                         "upper"=c(sapply(outcome.vars, function(i) boot.trajectory.eastern.placebo.cbw[[i]]$ci.upper),
#                                   sapply(outcome.vars, function(i) boot.trajectory.swiss.placebo.cbw[[i]]$ci.upper)),
#                         "Cluster"= c(rep(rep("Eastern",5), length(outcome.vars)), rep(rep("Swiss",5), length(outcome.vars))),
#                         "Outcome"=rep(outcome.vars,2*5))
# 
# # Plot
# cbw.placebo.plot <- ggplot(ci.values, aes(x=Outcome, y=t0, colour=as.factor(Cluster))) + 
#   geom_pointrange(aes(ymin=lower, ymax=upper))  +
#   labs(title="Retrospective prediction for later-treated", 
#        y=TeX("$S(\\hat{\\bar{\\tau}}_t)$"),
#        x="") + 
#   geom_hline(yintercept=0, linetype="dashed", color = "red") +
#   #ylim(-0.002, 0.002) +
# #  scale_x_discrete(labels=outcomes.labels, limits = levels(ci.values.m$Outcome)) +
#   # scale_shape_manual(name=expression(rho), values = c(1:5),
#   #                    labels = c(paste0(expression(rho),"=","1"),
#   #                               paste0(expression(rho),"=","2"),
#   #                               paste0(expression(rho),"=","3"),
#   #                               paste0(expression(rho),"=","4"),
#   #                               paste0(expression(rho),"=","5"))) +
#   scale_colour_manual(name="Cluster", values = c(  "Eastern" = wes_palette("Darjeeling1")[5],
#                                                    "Swiss" = wes_palette("Darjeeling1")[4]),
#                       labels=c("Eastern", "Swiss")) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme_set(theme_bw() + theme(legend.key=element_blank(), legend.title=element_text(size=10))) + theme(plot.title = element_text(hjust = 0.5, size=14)) + 
#   theme(axis.text.y = element_text(size=8)) + theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1, size=8))
# 
# ggsave(filename = paste0("plots/ci-cbw-",cf,".png"),plot = cbw.placebo.plot )