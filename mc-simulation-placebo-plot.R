############################################################################################
# Plot matrix completion simulation results                                                #
############################################################################################

library(tidyverse)
library(ggplot2)
library(data.table)
library(latex2exp)
library(dplyr)
library(grid)
library(gtable)
library(scales)

n.estimators <- 4

# Load results data

filenames <- c(list.files(path="outputs/20220123", pattern = ".rds", full.names = TRUE))

rmse.vars <- c('est_mc_weights_RMSE','est_model_ADH_RMSE','est_model_DID_RMSE','est_model_IFE_RMSE')
abs.bias.vars <- c('est_mc_weights_abs_bias','est_model_ADH_abs_bias','est_model_DID_abs_bias','est_model_IFE_abs_bias')
cp.vars <- c('est_mc_weights_cp','est_model_ADH_cp','est_model_DID_cp','est_model_IFE_cp')
boot.var.vars <- c('est_mc_weights_boot_var','est_model_ADH_boot_var','est_model_DID_boot_var','est_model_IFE_boot_var')

results <- list() # structure is: [[filename]][[metric]]
for(f in filenames){
  print(f)
  result.matrix <- readRDS(f)
  n <- nrow(result.matrix )
  rmse <- matrix(NA, n, n.estimators)
  abs.bias <- matrix(NA, n, n.estimators)
  CP <- matrix(NA, n, n.estimators)
  boot.var <- matrix(NA, n, n.estimators)
  for(i in rmse.vars){
    rmse[,which(rmse.vars==i)] <- unlist(result.matrix[,i])
  }
  for(i in abs.bias.vars){
    abs.bias[,which(abs.bias.vars==i)] <- unlist(result.matrix[,i])
  }
  for(i in cp.vars){
    CP[,which(cp.vars==i)] <- unlist(result.matrix[,i])
  }
  for(i in boot.var.vars){
    boot.var[,which(boot.var.vars==i)] <- unlist(result.matrix[,i])
  }
  results[[f]] <- list("rmse"=rmse,"abs_bias"=abs.bias,"CP"=CP,"boot_var"=boot.var,"n"=n)
}

# Create New lists
# structure is: [[estimator]][[filename]]

rmse <- list()
for(i in 1:length(rmse.vars)){
  rmse[[i]] <- lapply(1:length(filenames), function(f) results[[f]]$abs_bias[,i])
}

abs.bias <- list()
for(i in 1:length(abs.bias.vars)){
  abs.bias[[i]] <- lapply(1:length(filenames), function(f) results[[f]]$abs_bias[,i])
}

CP <- list()
for(i in 1:length(cp.vars)){
  CP[[i]] <- lapply(1:length(filenames), function(f) results[[f]]$CP[,i])
}

boot.var <- list()
for(i in 1:length(boot.var.vars)){
  boot.var[[i]] <- lapply(1:length(filenames), function(f) results[[f]]$boot_var[,i])
}

# Create dataframe for plot
results.df <- data.frame("rmse"=as.numeric(unlist(rmse)),
                        "abs_bias"=as.numeric(unlist(abs.bias)),
                         "Coverage"=as.numeric(unlist(CP)),
                         "boot_var"=as.numeric(unlist(boot.var)),
                         "Estimator"=c(rep("MC (weights)",length.out=length(c(unlist(CP[[1]])))),
                         rep("SCM",length.out=length(c(unlist(CP[[2]])))),
                         rep("DID",length.out=length(c(unlist(CP[[3]])))),
                         rep("IFE",length.out=length(c(unlist(CP[[4]]))))),
                         "filename"=c(rep(unlist(sapply(1:length(filenames), function(i) rep(filenames[i], length.out=n))), n.estimators)))

results.df$N_t <- NA
results.df$T0 <- NA
 
for(s in c("T0_0.5","T0_0.65","T0_0.8","T0_0.95")){
  if(length(grep(s, results.df$filename))>0){
    print(s)
    if(s=="T0_0.5"){
      T0 <- "0.5"
    }
    if(s=="T0_0.65"){
      T0 <- "0.65"
    }
    if(s=="T0_0.80"){
      T0 <- "0.80"
    }
    if(s=="T0_0.95"){
      T0 <- "0.95"
    }
    results.df[grep(s, results.df$filename),]$T0 <- T0
  }
}

for(s in c("N_t0.5")){
  if(length(grep(s, results.df$filename))>0){
    print(s)
    if(s=="N_t0.5"){
      N_t <- "0.5"
    }
    results.df[grep(s, results.df$filename),]$N_t  <- N_t
  }
}

# create coverage rate variable

results.df <- results.df %>% 
  group_by(Estimator,T0,N_t) %>% 
  mutate(CP = mean(Coverage)) 

# reshape and plot
results.df$id <- with(results.df, paste(T0,N_t, sep = "_"))
results_long <- reshape2::melt(results.df[!colnames(results.df) %in% c("id","filename")], id.vars=c("Estimator","T0","N_t"))  # convert to long format

# rmse (NxT)
sim.results.rmse <- ggplot(data=results_long[results_long$variable=="rmse",],
                               aes(x=factor(T0), y=value, fill=Estimator))  + geom_boxplot(outlier.alpha = 0.3,outlier.size = 1, outlier.stroke = 0.1, lwd=0.25) +
 ylab("RMSE") +  xlab(TeX('Placebo $(a_i^{\\prime}/T)$')) +
  scale_fill_discrete(name = "Estimator:") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA),breaks= pretty_breaks())+
  theme(legend.position="right") +   theme(plot.title = element_text(hjust = 0.5, family="serif", size=16)) +
  theme(axis.title=element_text(family="serif", size=16)) +
  theme(axis.text.y=element_text(family="serif", size=14)) +
  theme(axis.text.x=element_text(family="serif", size=14)) +
  theme(legend.text=element_text(family="serif", size = 14)) +
  theme(legend.title=element_text(family="serif", size = 14)) +
  theme(strip.text.x = element_text(family="serif", size = 14)) +
  theme(strip.text.y = element_text(family="serif", size = 14)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l =0))) +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l =0))) +
  theme(panel.spacing = unit(1, "lines")) 

ggsave("plots/simulation_rmse.png",plot = sim.results.rmse)

# abs.bias (NxT)
sim.results.abs.bias <- ggplot(data=results_long[results_long$variable=="abs_bias",],
                           aes(x=factor(T0), y=value, fill=Estimator))  + geom_boxplot(outlier.alpha = 0.3,outlier.size = 1, outlier.stroke = 0.1, lwd=0.25) +
 ylab("Absolute bias") +  xlab(TeX('Placebo $(a_i^{\\prime}/T)$')) +
  scale_fill_discrete(name = "Estimator:") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA),breaks= pretty_breaks())+
  theme(legend.position="right") +   theme(plot.title = element_text(hjust = 0.5, family="serif", size=16)) +
  theme(axis.title=element_text(family="serif", size=16)) +
  theme(axis.text.y=element_text(family="serif", size=14)) +
  theme(axis.text.x=element_text(family="serif", size=14)) +
  theme(legend.text=element_text(family="serif", size = 14)) +
  theme(legend.title=element_text(family="serif", size = 14)) +
  theme(strip.text.x = element_text(family="serif", size = 14)) +
  theme(strip.text.y = element_text(family="serif", size = 14)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l =0))) +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l =0))) +
  theme(panel.spacing = unit(1, "lines")) 

ggsave("plots/simulation_abs_bias.png",plot = sim.results.abs.bias)

# coverage
sim.results.coverage <- ggplot(data=results_long[results_long$variable=="CP",],
                           aes(x=factor(T0), y=value, colour=Estimator, group=forcats::fct_rev(Estimator)))  +   geom_line()  +
  xlab(TeX('Placebo $(a_i^{\\prime}/T)$')) + ylab("Coverage probability (%)") +
  scale_colour_discrete(name = "Estimator:") +
  scale_y_continuous(breaks= pretty_breaks()) +
  geom_hline(yintercept = 0.95, linetype="dotted")+
  theme(legend.position="right") +   theme(plot.title = element_text(hjust = 0.5, family="serif", size=16)) +
  theme(axis.title=element_text(family="serif", size=16)) +
  theme(axis.text.y=element_text(family="serif", size=14)) +
  theme(axis.text.x=element_text(family="serif", size=14)) +
  theme(legend.text=element_text(family="serif", size = 14)) +
  theme(legend.title=element_text(family="serif", size = 14)) +
  theme(strip.text.x = element_text(family="serif", size = 14)) +
  theme(strip.text.y = element_text(family="serif", size = 14)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l =0))) +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l =0))) +
  theme(panel.spacing = unit(1, "lines"))

ggsave("plots/simulation_coverage.png",plot = sim.results.coverage)

# boot_var

sim.results.boot.var <- ggplot(data=results_long[results_long$variable=="boot_var",],
                           aes(x=factor(T0), y=value, fill=Estimator))  + geom_boxplot(outlier.alpha = 0.3,outlier.size = 1, outlier.stroke = 0.1, lwd=0.25) +
  xlab(TeX('Placebo $(a_i^{\\prime}/T)$'))  + ylab("Bootstrap variance") +
  scale_fill_discrete(name = "Estimator:") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA),breaks= pretty_breaks())+
  theme(legend.position="right") +   theme(plot.title = element_text(hjust = 0.5, family="serif", size=16)) +
  theme(axis.title=element_text(family="serif", size=16)) +
  theme(axis.text.y=element_text(family="serif", size=14)) +
  theme(axis.text.x=element_text(family="serif", size=14)) +
  theme(legend.text=element_text(family="serif", size = 14)) +
  theme(legend.title=element_text(family="serif", size = 14)) +
  theme(strip.text.x = element_text(family="serif", size = 14)) +
  theme(strip.text.y = element_text(family="serif", size = 14)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l =0))) +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l =0))) +
  theme(panel.spacing = unit(1, "lines"))

ggsave("plots/simulation_boot_var.png",plot = sim.results.boot.var)

# Get color hues

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_color_hue(n.estimators)