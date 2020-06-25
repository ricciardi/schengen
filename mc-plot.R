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

source("TsPlot.R")

PlotMCCapacity <- function(observed,main,y.title,mc_est,boot_result,treated,control,vline,breaks,labels,att.label,rev){
  ## Create time series data
  
  predicted <- mc_est$Mhat
  
  if(rev){
    pointwise <- abs(mc_est$impact)
  }
  else{
    pointwise <- mc_est$impact
  }
  
  pointwise.se <- matrix(apply(boot_result$t, 2, sd), nrow=dim(pointwise)[1], ncol=dim(pointwise)[2], byrow=FALSE)
  
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
  
  ts.plot <- TsPlot(df=ts.means.m,main=main, y.title=y.title,vline,breaks,labels,rev)
  
  return(ts.plot)
}

## Plot time-series

outcome.vars <- c("CBWbord","CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_0","seekdur_1_2","seekdur_3more")
outcomes.labels <- c("Share of residents working in border region",
                     "Share of employed residents working in border region",
                     "Regional employment rate",
                     "Average total working hours",
                     "Unemployment rate",
                     "Inactivity rate",
                     "% unemployed for < 1 month",
                     "% unemployed for < 1-2 months",
                     "% unemployed for < 3 months")

covarflag <- c("","-covars")

for(o in outcome.vars){
  for(c in covarflag){
    print(c)
    print(o)
    
    ## Analysis 1: ST vs AT (retrospective, X=CBW) 
    
    outcomes.cbw <- readRDS(paste0("data/outcomes-cbw",o,".rds"))
    mc.estimates.cbw <- readRDS(paste0("results/mc-estimates-cbw",o,c,".rds"))
    boot.cbw <- readRDS(paste0("results/boot-cbw",o,c,".rds"))
    
    mc.plot <- PlotMCCapacity(observed = outcomes.cbw$M, 
                              y.title=outcomes.labels[which(outcome.vars==o)],
                              main = "Retrospective prediction for later-treated in Eastern cluster",
                              mc_est=mc.estimates.cbw, boot_result=boot.cbw, 
                              treated=outcomes.cbw$treated, control=outcomes.cbw$control, vline=20091,vline2=20111,
                              breaks=c(20042,20091,20111,20121,20161,20184),
                              labels=c(20042,20091,20111,20121,20161,20184),
                              att.label = TeX("$\\hat{\\bar{\\tau}}$"),
                              rev=TRUE)
    
    ggsave(paste0("plots/mc-estimates-cbw",o,c,".png"), mc.plot, width=8.5, height=11)
    
    ## Analysis 2:  ST vs NT (forward, X=LM)
    
    outcomes.lm <- readRDS(paste0("data/outcomes-lm",o,".rds"))
    mc.estimates.lm <- readRDS(paste0("results/mc-estimates-lm",o,c,".rds"))
    boot.lm <- readRDS(paste0("results/boot-lm",o,c,".rds"))
    
    mc.plot <- PlotMCCapacity(observed = outcomes.lm$M, 
                              y.title=outcomes.labels[which(outcome.vars==o)],
                              main = "Prospective prediction for later-treated in Eastern cluster",
                              mc_est=mc.estimates.lm, boot_result=boot.lm, 
                              treated=outcomes.lm$treated, control=outcomes.lm$control, vline=20091,vline2=20111,
                              breaks=c(20042,20091,20111,20121,20161,20184),
                              labels=c(20042,20091,20111,20121,20161,20184),
                              att.label = TeX("$\\hat{\\bar{\\tau}}$"),
                              rev=FALSE)
    
    ggsave(paste0("plots/mc-estimates-lm",o,c,".png"), mc.plot, width=8.5, height=11)
    
  }
}

## Plot p-values

for(c in covarflag){
  print(c)
  ## Analysis 1: ST vs AT (retrospective, X=CBW) 
  
  iid <- lapply(outcome.vars, function(o){
    p <- readRDS(paste0("results/iid-cbw",o,c,".rds"))
    return(p)
  })
  names(iid) <- outcome.vars
  
  iid.block <- lapply(outcome.vars, function(o){
    p <- readRDS(paste0("results/iid-block-cbw",o,c,".rds"))
    return(p)
  })
  names(iid.block) <- outcome.vars
  
  moving.block <- lapply(outcome.vars, function(o){
    p <- readRDS(paste0("results/moving-block-cbw",o,c,".rds"))
    return(p)
  })
  names(moving.block) <- outcome.vars
  
  p.values <- data.frame("iid"=c(sapply(outcome.vars, function(i) iid[[i]]$p)[1,],sapply(outcome.vars, function(i) iid[[i]]$p)[2,]),
                         "iid.block"=c(sapply(outcome.vars, function(i) iid.block[[i]]$p)[1,],sapply(outcome.vars, function(i) iid.block[[i]]$p)[2,]),
                         "moving.block"=c(sapply(outcome.vars, function(i) moving.block[[i]]$p)[1,],sapply(outcome.vars, function(i) moving.block[[i]]$p)[2,]),
                         "q"=c(rep(1, length(outcome.vars)), rep(2, length(outcome.vars))),
                         "outcomes"=c(rep(outcome.vars, length(outcome.vars)), rep(outcome.vars, length(outcome.vars))))
  
  p.values.m <-melt(p.values,id.vars = c("outcomes","q"))
  
  p.values.m$Type <- paste0(p.values.m$variable, ", q=", p.values.m$q)
  
  # Plot
  mc.pvals.plot <- ggplot(p.values.m, aes(x=outcomes, y=value)) + 
    geom_point(stat='identity', aes(col=variable,shape=factor(q)), size=3, alpha=0.5)  +
    labs(title="Retrospective prediction for later-treated", 
         y="Randomization p-values",
         y="Outcomes") + 
    geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
    scale_y_continuous(breaks=c(0.05,0.25,0.5,0.75,1), 
                       labels=c("0.05","0.25","0.50","0.75","1")) +
    scale_x_discrete(labels=rev(outcomes.labels), limits = rev(levels(p.values.m$outcomes)))+
    coord_flip() +
    scale_shape_manual(name="", values = c("1" = 2,
                                           "2" = 4),
                       labels=c("q=1", "q=2")) +
    scale_colour_manual(name="Randomization type", values = c(  "iid" = wes_palette("Darjeeling1")[1],
                                                                "iid.block" = wes_palette("Darjeeling1")[2], 
                                                                "moving.block" = wes_palette("Darjeeling1")[4]),
                        labels=c("IID", "IID Block", 
                                 "Moving block")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme_set(theme_bw() + theme(legend.key=element_blank(), legend.title=element_text(size=10))) + theme(plot.title = element_text(hjust = 0.5, size=14)) + 
    theme(axis.text.y = element_text(size=8))
  
  ggsave(filename = "plots/pvals-cbw.png",plot = mc.pvals.plot)
  
  ## Analysis 2:  ST vs NT (forward, X=LM)
  
  iid <- lapply(outcome.vars, function(o){
    p <- readRDS(paste0("results/iid-lm",o,c,".rds"))
    return(p)
  })
  names(iid) <- outcome.vars
  
  iid.block <- lapply(outcome.vars, function(o){
    p <- readRDS(paste0("results/iid-block-lm",o,c,".rds"))
    return(p)
  })
  names(iid.block) <- outcome.vars
  
  moving.block <- lapply(outcome.vars, function(o){
    p <- readRDS(paste0("results/moving-block-lm",o,c,".rds"))
    return(p)
  })
  names(moving.block) <- outcome.vars
  
  p.values <- data.frame("iid"=c(sapply(outcome.vars, function(i) iid[[i]]$p)[1,],sapply(outcome.vars, function(i) iid[[i]]$p)[2,]),
                         "iid.block"=c(sapply(outcome.vars, function(i) iid.block[[i]]$p)[1,],sapply(outcome.vars, function(i) iid.block[[i]]$p)[2,]),
                         "moving.block"=c(sapply(outcome.vars, function(i) moving.block[[i]]$p)[1,],sapply(outcome.vars, function(i) moving.block[[i]]$p)[2,]),
                         "q"=c(rep(1, length(outcome.vars)), rep(2, length(outcome.vars))),
                         "outcomes"=c(rep(outcome.vars, length(outcome.vars)), rep(outcome.vars, length(outcome.vars))))
  
  p.values.m <-melt(p.values,id.vars = c("outcomes","q"))
  
  p.values.m$Type <- paste0(p.values.m$variable, ", q=", p.values.m$q)
  
  # Plot
  mc.pvals.plot <- ggplot(p.values.m, aes(x=outcomes, y=value)) + 
    geom_point(stat='identity', aes(col=variable,shape=factor(q)), size=3, alpha=0.5)  +
    labs(title="Prospective prediction for later-treated", 
         y="Randomization p-values",
         x="Outcomes") + 
    geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
    scale_y_continuous(breaks=c(0.05,0.25,0.5,0.75,1), 
                       labels=c("0.05","0.25","0.50","0.75","1")) +
    scale_x_discrete(labels=rev(outcomes.labels), limits = rev(levels(p.values.m$outcomes)))+
    coord_flip() +
    scale_shape_manual(name="", values = c("1" = 2,
                                           "2" = 4),
                       labels=c("q=1", "q=2")) +
    scale_colour_manual(name="Randomization type", values = c(  "iid" = wes_palette("Darjeeling1")[1],
                                                                "iid.block" = wes_palette("Darjeeling1")[2], 
                                                                "moving.block" = wes_palette("Darjeeling1")[4]),
                        labels=c("IID", "IID Block", 
                                 "Moving block")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme_set(theme_bw() + theme(legend.key=element_blank(), legend.title=element_text(size=10))) + theme(plot.title = element_text(hjust = 0.5, size=14)) + 
    theme(axis.text.y = element_text(size=8))
  
  ggsave(filename = "plots/pvals-lm.png",plot = mc.pvals.plot)
}