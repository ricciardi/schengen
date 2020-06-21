###################################
# Plot Observed/predicted and causal impacts #
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

PlotMCCapacity <- function(observed,main,y.title,t0,mc_est,boot_result,treated,control,vline,breaks,labels,att.label,rev){
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
  ts.means.m$series[grep("observed.", ts.means.m$variable)] <- "Observed/predicted"
  ts.means.m$series[grep("predicted.", ts.means.m$variable)] <- "Observed/predicted"
  ts.means.m$series[grep("pointwise.", ts.means.m$variable)] <- att.label
  
  ts.means.m$series<- factor(ts.means.m$series, levels=c("Observed/predicted", att.label)) # reverse order
  
  ts.plot <- TsPlot(df=ts.means.m,main=main, y.title=y.title,vline,breaks,labels,rev)
  
  return(ts.plot)
}

## Plot time-series

outcome.vars <- c("CBWbord","CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_0","seekdur_1_2","seekdur_3more")
outcomes.labels <- c("Share of residents working in a border region",
                     "Share of residents working in a border region, conditional on employment",
                     "Regional employment rate",
                     "Average total working hours",
                     "Unemployment rate",
                     "Inactivity rate",
                     "% of unemployed with unemployment duration less than 1 month",
                     "% of unemployed with unemployment duration less than 1-2 months",
                     "% of unemployed with unemployment duration less than 3 months or more")

for(o in outcome.vars){
  print(o)
  
  ## Analysis 1: ST vs AT (retrospective, X=CBW) 
  
  # Eastern cluster
  
  outcomes.cbw.eastern <- readRDS(paste0("data/outcomes-cbw-eastern-",o,".rds"))
  mc.estimates.cbw.eastern <- readRDS(paste0("results/mc-estimates-cbw-eastern-",o,".rds"))
  boot.cbw.eastern <- readRDS(paste0("results/boot-cbw-eastern-",o,".rds"))
  
  mc.plot <- PlotMCCapacity(observed = outcomes.cbw.eastern$M, 
                            y.title=outcomes.labels[which(outcome.vars==o)],
                            main = "Retrospective prediction for later-treated",
                            t0=which(colnames(outcomes.cbw.eastern$M)=="20081"),
                            mc_est=mc.estimates.cbw.eastern, boot_result=boot.cbw.eastern, 
                            treated=outcomes.cbw.eastern$treated, control=outcomes.cbw.eastern$control, vline=20081,
                            breaks=c(20042,20081,20121,20161,20184),
                            labels=c(20042,20081,20121,20161,20184),
                            att.label = "ATT",
                            rev=TRUE)
  
  ggsave(paste0("plots/mc-estimates-cbw-eastern-",o,".png"), mc.plot, width=8.5, height=11)
  
  # Swiss cluster
  
  outcomes.cbw.swiss <- readRDS(paste0("data/outcomes-cbw-swiss-",o,".rds"))
  mc.estimates.cbw.swiss <- readRDS(paste0("results/mc-estimates-cbw-swiss-",o,".rds"))
  boot.cbw.swiss <- readRDS(paste0("results/boot-cbw-swiss-",o,".rds"))
  
  mc.plot <- PlotMCCapacity(observed = outcomes.cbw.swiss$M, 
                            y.title=outcomes.labels[which(outcome.vars==o)],
                            main = "Retrospective prediction for later-treated",
                            t0=which(colnames(outcomes.cbw.swiss$M)=="20072"),
                            mc_est=mc.estimates.cbw.swiss, boot_result=boot.cbw.swiss, 
                            treated=outcomes.cbw.swiss$treated, control=outcomes.cbw.swiss$control, vline=20072,
                            breaks=c(20042,20081,20121,20161,20184),
                            labels=c(20042,20081,20121,20161,20184),
                            att.label = "ATT",
                            rev=TRUE)
  
  ggsave(paste0("plots/mc-estimates-cbw-swiss-",o,".png"), mc.plot, width=8.5, height=11)
  
  ## Analysis 2:  ST vs NT (forward, X=LM)
  
  # Eastern cluster
  
  outcomes.lm.eastern <- readRDS(paste0("data/outcomes-lm-eastern-",o,".rds"))
  mc.estimates.lm.eastern <- readRDS(paste0("results/mc-estimates-lm-eastern-",o,".rds"))
  boot.lm.eastern <- readRDS(paste0("results/boot-lm-eastern-",o,".rds"))
  
  mc.plot <- PlotMCCapacity(observed = outcomes.lm.eastern$M, 
                            y.title=outcomes.labels[which(outcome.vars==o)],
                            main = "Prospective prediction for later-treated",
                            t0=which(colnames(outcomes.lm.eastern$M)=="20081"),
                            mc_est=mc.estimates.lm.eastern, boot_result=boot.lm.eastern, 
                            treated=outcomes.lm.eastern$treated, control=outcomes.lm.eastern$control, vline=20081,
                            breaks=c(20042,20081,20121,20161,20184),
                            labels=c(20042,20081,20121,20161,20184),
                            att.label = "ATT",
                            rev=FALSE)
  
  ggsave(paste0("plots/mc-estimates-lm-eastern-",o,".png"), mc.plot, width=8.5, height=11)
  
  # Swiss cluster
  
  outcomes.lm.swiss <- readRDS(paste0("data/outcomes-lm-swiss-",o,".rds"))
  mc.estimates.lm.swiss <- readRDS(paste0("results/mc-estimates-lm-swiss-",o,".rds"))
  boot.lm.swiss <- readRDS(paste0("results/boot-lm-swiss-",o,".rds"))
  
  mc.plot <- PlotMCCapacity(observed = outcomes.lm.swiss$M, 
                            y.title=outcomes.labels[which(outcome.vars==o)],
                            main = "Prospective prediction for later-treated",
                            t0=which(colnames(outcomes.lm.swiss$M)=="20072"),
                            mc_est=mc.estimates.lm.swiss, boot_result=boot.lm.swiss, 
                            treated=outcomes.lm.swiss$treated, control=outcomes.lm.swiss$control, vline=20072,
                            breaks=c(20042,20081,20121,20161,20184),
                            labels=c(20042,20081,20121,20161,20184),
                            att.label = "ATT",
                            rev=FALSE)
  
  ggsave(paste0("plots/mc-estimates-lm-swiss-",o,".png"), mc.plot, width=8.5, height=11)
}

## Plot p-values ## REVISE

iid.placebo <- readRDS("results/iid_placebo.rds")

iid.block.placebo <- readRDS("results/iid_block_placebo.rds")

moving.block.placebo <- readRDS("results/moving_block_placebo.rds")

p.values <- data.frame("iid"=c(sapply(taus, function(i) iid.placebo[[i]]$p)[1,],sapply(taus, function(i) iid.placebo[[i]]$p)[2,]),
                       "iid.block"=c(sapply(taus, function(i) iid.block.placebo[[i]]$p)[1,],sapply(taus, function(i) iid.block.placebo[[i]]$p)[2,]),
                       "moving.block"=c(sapply(taus, function(i) moving.block.placebo[[i]]$p)[1,],sapply(taus, function(i) moving.block.placebo[[i]]$p)[2,]),
                       "q"=c(rep(1, length(taus)), rep(2, length(taus))),
                       "tau"=rep(taus,2))

p.values.m <-melt(p.values,id.vars = c("tau","q"))

p.values.m$Type <- paste0(p.values.m$variable, ", q=", p.values.m$q)

# Plot
mc.placebo.plot <- ggplot(p.values.m, aes(x=value, y=tau)) + 
  geom_point(stat='identity', aes(col=variable,shape=factor(q)), size=3, alpha=0.5)  +
  labs(title="Placebo test", 
       x="Randomization p-values",
       y=expression(tau)) + 
  geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
  scale_x_continuous(breaks=c(0,0.05,0.1,0.25,0.5,0.75,1), 
                     labels=c("0","0.05","0.10","0.25","0.50","0.75","1")) +
  scale_y_continuous(breaks=seq(3,15,3), 
                     labels=c("3",  "6", "9", "12","15"))+
  coord_flip() +
  scale_shape_manual(name="", values = c("1" = 2,
                                         "2" = 4),
                     labels=c("q=1", "q=2")) +
  scale_colour_manual(name="Randomization methods", values = c(  "iid" = wes_palette("Darjeeling1")[1],
                                                                 "iid.block" = wes_palette("Darjeeling1")[2], 
                                                                 "moving.block" = wes_palette("Darjeeling1")[4]),
                      labels=c("IID", "IID Block", 
                               "Moving block")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme_set(theme_bw() + theme(legend.key=element_blank())) + theme(plot.title = element_text(hjust = 0.5)) # rm background

ggsave(filename = "plots/mc-placebo-plot.png",plot = mc.placebo.plot)