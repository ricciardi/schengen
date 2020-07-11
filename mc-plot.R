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
  
  ts.means.m <- ts.means.m %>% # pivot CI
    mutate(upper = 2*value -boot.lower,
           lower = 2*value -boot.upper)
  
  # Labels
  
  ts.means.m$series <- NA
  ts.means.m$series[grep("observed.", ts.means.m$variable)] <- "Time-series"
  ts.means.m$series[grep("predicted.", ts.means.m$variable)] <- "Time-series"
  ts.means.m$series[grep("pointwise.", ts.means.m$variable)] <- att.label
  
  ts.means.m$series<- factor(ts.means.m$series, levels=c("Time-series", att.label)) # reverse order
  
  ts.plot <- TsPlot(df=ts.means.m,main=main, y.title=y.title,vline,vline2,breaks,labels,rev)
  
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
    
  #  if(cf=="" && o %in% c("N_CBWbord")) next

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
  }
} 

# ## Plot p-values
# 
# outcome.vars <- c("CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_0","seekdur_1_2","seekdur_3more")
# for(o in outcome.vars){
#   for(c in covarflag){
#     print(c)
#     
#     iid.block.eastern <- lapply(outcome.vars, function(o){
#       p <- readRDS(paste0("results/iid-block-cbw-eastern-",o,cf,".rds"))
#       return(p)
#     })
#     
#     names(iid.block.eastern) <- outcome.vars
#     
#     iid.block.swiss<- lapply(outcome.vars, function(o){
#       p <- readRDS(paste0("results/iid-block-cbw-swiss-",o,cf,".rds"))
#       return(p)
#     })
#     
#     names(iid.block.swiss) <- outcome.vars
#     
#     p.values <- data.frame("p"=c(sapply(outcome.vars, function(i) iid.block.eastern[[i]]$p), sapply(outcome.vars, function(i) iid.block.swiss[[i]]$p)),
#                            "clusters"=c(rep("Eastern", each=length(outcome.vars)), rep("Swiss", each=length(outcome.vars))),
#                            "outcomes"=c(rep(outcome.vars, times=2)))
#     
#     p.values.m <-melt(p.values,id.vars = c("outcomes","clusters"))
#     
# 
#     # Plot
#     mc.pvals.plot <- ggplot(p.values.m, aes(x=outcomes, y=value)) + 
#       geom_point(stat='identity', aes(col=variable,shape=factor(q)), size=3, alpha=0.5)  +
#       labs(title="Retrospective prediction for later-treated", 
#            y="Randomization p-values",
#            y="Outcomes") + 
#       geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
#       scale_y_continuous(breaks=c(0.05,0.25,0.5,0.75,1), 
#                          labels=c("0.05","0.25","0.50","0.75","1")) +
#       scale_x_discrete(labels=rev(outcomes.labels), limits = rev(levels(p.values.m$outcomes)))+
#       coord_flip() +
#       scale_shape_manual(name="", values = c("1" = 2,
#                                              "2" = 4),
#                          labels=c("q=1", "q=2")) +
#       scale_colour_manual(name="Randomization type", values = c(  "iid" = wes_palette("Darjeeling1")[1],
#                                                                   "iid.block" = wes_palette("Darjeeling1")[2], 
#                                                                   "moving.block" = wes_palette("Darjeeling1")[4]),
#                           labels=c("IID", "IID Block", 
#                                    "Moving block")) +
#       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#             panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme_set(theme_bw() + theme(legend.key=element_blank(), legend.title=element_text(size=10))) + theme(plot.title = element_text(hjust = 0.5, size=14)) + 
#       theme(axis.text.y = element_text(size=8))
#     
#     ggsave(filename = "plots/pvals-cbw-.png",plot = mc.pvals.plot)
#   }
# }