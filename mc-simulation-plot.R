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

n.estimators <- 7

# Load results data
output_dir <- "mc-simulations"

filenames <- c(list.files(path="outputs/20220103", pattern = ".rds", full.names = TRUE),
               list.files(path="outputs/20211229", pattern = ".rds", full.names = TRUE))

abs.bias.vars <- c('est_mc_plain_abs_bias','est_mc_weights_abs_bias','est_mc_covars_abs_bias','est_mc_weights_covars_abs_bias','est_model_ADH_abs_bias','est_model_DID_abs_bias','est_model_IFE_abs_bias')
rel.abs.bias.vars <- c('est_mc_plain_rel_abs_bias','est_mc_weights_rel_abs_bias','est_mc_covars_rel_abs_bias','est_mc_weights_covars_rel_abs_bias','est_model_ADH_rel_abs_bias','est_model_DID_rel_abs_bias','est_model_IFE_rel_abs_bias')
cp.vars <- c('est_mc_plain_cp','est_mc_weights_cp','est_mc_covars_cp','est_mc_weights_covars_cp','est_model_ADH_cp','est_model_DID_cp','est_model_IFE_cp')
boot.var.vars <- c('est_mc_plain_boot_var','est_mc_weights_boot_var','est_mc_covars_boot_var','est_mc_weights_covars_boot_var','est_model_ADH_boot_var','est_model_DID_boot_var','est_model_IFE_boot_var')

results <- list() # structure is: [[filename]][[metric]]
for(f in filenames){
  print(f)
  result.matrix <- readRDS(f)
  n <- nrow(result.matrix )
  abs.bias <- matrix(NA, n, n.estimators)
  rel.abs.bias <- matrix(NA, n, n.estimators)
  CP <- matrix(NA, n, n.estimators)
  boot.var <- matrix(NA, n, n.estimators)
  for(i in abs.bias.vars){
    abs.bias[,which(abs.bias.vars==i)] <- unlist(result.matrix[,i])
  }
  for(i in rel.abs.bias.vars){
    rel.abs.bias[,which(rel.abs.bias.vars==i)] <- unlist(result.matrix[,i])
  }
  for(i in cp.vars){
    CP[,which(cp.vars==i)] <- unlist(result.matrix[,i])
  }
  for(i in boot.var.vars){
    boot.var[,which(boot.var.vars==i)] <- unlist(result.matrix[,i])
  }
  fr_obs <- unlist(result.matrix[,"fr_obs"])
  results[[f]] <- list("abs_bias"=abs.bias,"rel_abs_bias"=rel.abs.bias,"CP"=CP,"boot_var"=boot.var,"fr_obs"=fr_obs,"n"=n)
}

# Create New lists
# structure is: [[estimator]][[filename]]

abs.bias <- list()
for(i in 1:length(abs.bias.vars)){
  abs.bias[[i]] <- lapply(1:length(filenames), function(f) results[[f]]$abs_bias[,i])
}

rel.abs.bias <- list()
for(i in 1:length(rel.abs.bias.vars)){
  rel.abs.bias[[i]] <- lapply(1:length(filenames), function(f) results[[f]]$rel_abs_bias[,i])
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
results.df <- data.frame("abs_bias"=as.numeric(unlist(abs.bias)),
                         "rel_abs_bias"=as.numeric(unlist(rel.abs.bias)),
                         "Coverage"=as.numeric(unlist(CP)),
                         "boot_var"=as.numeric(unlist(boot.var)),
                         "fr_obs"= rep(sapply(1:length(filenames), function(f) results[[f]]$fr_obs), n.estimators),
                         "Estimator"=c(rep("MC (Athey et al.)",length.out=length(c(unlist(CP[[1]])))), 
                                       rep("MC (+ weights)",length.out=length(c(unlist(CP[[2]])))),
                                       rep("MC (+ covars.)",length.out=length(c(unlist(CP[[2]])))),
                         rep("MC (weights + covars.)",length.out=length(c(unlist(CP[[3]])))),
                         rep("SCM",length.out=length(c(unlist(CP[[4]])))),
                         rep("DID",length.out=length(c(unlist(CP[[5]])))),
                         rep("IFE",length.out=length(c(unlist(CP[[5]]))))),
                         "filename"=c(rep(unlist(sapply(1:length(filenames), function(i) rep(filenames[i], length.out=n))), n.estimators)))

results.df$NT <- NA
results.df$R <- NA
results.df$shift_sc <- NA
results.df$T0 <- NA
  
for(s in c("N_40_T_40","N_50_T_50","N_60_T_60")){
  if(length(grep(s, results.df$filename))>0){
    print(s)
    if(s=="N_40_T_40"){
      NT <- "1600"
    }
    if(s=="N_50_T_50"){
      NT <- "2500"
    }
    if(s=="N_60_T_60"){
      NT <- "3600"
    }
    results.df[grep(s, results.df$filename),]$NT <- NT
  }
}

for(s in c("R_20","R_25","R_30")){
  if(length(grep(s, results.df$filename))>0){
    print(s)
    if(s=="R_20"){
      R <- "10"
    }
    if(s=="R_25"){
      R <- "20"
    }
    if(s=="R_30"){
      R <- "40"
    }
    results.df[grep(s, results.df$filename),]$R <- R
  }
}

for(s in c("T0_20","T0_26","T0_32","T0_39")){
  if(length(grep(s, results.df$filename))>0){
    print(s)
    if(s=="T0_20"){
      T0 <- "0.5"
    }
    if(s=="T0_26"){
      T0 <- "0.65"
    }
    if(s=="T0_32"){
      T0 <- "0.8"
    }
    if(s=="T0_39"){
      T0 <- "0.975"
    }
    results.df[grep(s, results.df$filename),]$T0 <- T0
  }
}

for(s in c("shift_sc_0.1")){
  if(length(grep(s, results.df$filename))>0){
    print(s)
    if(s=="shift_sc_0.1"){
      shift_sc <- "0.1"
    }
    results.df[grep(s, results.df$filename),]$shift_sc  <- shift_sc
  }
}

# create coverage rate variable

results.df <- results.df %>% 
  group_by(Estimator,NT,R,T0,shift_sc) %>% 
  mutate(CP = mean(Coverage)) 

# reshape and plot
results.df$id <- with(results.df, paste(NT,R,T0,shift_sc, sep = "_"))
results_long <- reshape2::melt(results.df[!colnames(results.df) %in% c("id","filename")], id.vars=c("Estimator","NT","R","T0","shift_sc","fr_obs"))  # convert to long format

variable_names <- list(
  '1600'= TeX("$40 \\times 40"),
  '3600'= TeX("$50 \\times 50"),
  '6400'= TeX("$60 \\times 60")) 

labeller <- function(variable,value){
  return(variable_names[value])
}

# abs.bias (NxT)
sim.results.abs.bias <- ggplot(data=results_long[results_long$variable=="abs_bias",],
                           aes(x=factor(T0), y=value, fill=Estimator))  + geom_boxplot(outlier.alpha = 0.3,outlier.size = 1, outlier.stroke = 0.1, lwd=0.25) +
  facet_grid(shift_sc ~  NT, scales = "free", labeller=labeller)  + ylab("Absolute bias") +  xlab(TeX('Placebo $(a_i^{\\prime}/T)$')) +
  scale_fill_discrete(name = "Estimator:") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(breaks= pretty_breaks()) +
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

# Get the ggplot grob
z.abs.bias <- ggplotGrob(sim.results.abs.bias)

# Labels 
labelR <- " "
labelT <- TeX("$N \\times T")

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(z.abs.bias$layout, grepl("strip-r", name), select = t:r)
posT <- subset(z.abs.bias$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- z.abs.bias$widths[max(posR$r)]    # width of current right strips
height <- z.abs.bias$heights[min(posT$t)]  # height of current top strips

z.abs.bias <- gtable_add_cols(z.abs.bias, width, max(posR$r))  
z.abs.bias <- gtable_add_rows(z.abs.bias, height, min(posT$t)-1)

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, fill = "grey85")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 16, col = "grey10"))))

stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, fill = "grey85")),
  textGrob(labelT, gp = gpar(fontsize = 16, col = "grey10"))))

# Position the grobs in the gtable
z.abs.bias <- gtable_add_grob(z.abs.bias, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
z.abs.bias <- gtable_add_grob(z.abs.bias, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
z.abs.bias <- gtable_add_cols(z.abs.bias, unit(1/5, "line"), max(posR$r))
z.abs.bias <- gtable_add_rows(z.abs.bias, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(z.abs.bias)

ggsave("plots/simulation_abs_bias.png",plot = z.abs.bias,scale=2)

# rel.abs.bias (NxT)
sim.results.rel.abs.bias <- ggplot(data=results_long[results_long$variable=="rel_abs_bias",],
                               aes(x=factor(T0), y=value, fill=Estimator))  + geom_boxplot(outlier.alpha = 0.3,outlier.size = 1, outlier.stroke = 0.1, lwd=0.25) +
  facet_grid(shift_sc ~  NT, scales = "free", labeller=labeller)  + ylab("Relative absolute bias") +  xlab(TeX('Placebo $(a_i^{\\prime}/T)$')) +
  scale_fill_discrete(name = "Estimator:") +
  scale_y_continuous(breaks= pretty_breaks()) +
  coord_cartesian(ylim = c(0,5)) +
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

# Get the ggplot grob
z.rel.abs.bias <- ggplotGrob(sim.results.rel.abs.bias)

# Labels 
labelR <- " "
labelT <- TeX("$N \\times T")

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(z.rel.abs.bias$layout, grepl("strip-r", name), select = t:r)
posT <- subset(z.rel.abs.bias$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- z.rel.abs.bias$widths[max(posR$r)]    # width of current right strips
height <- z.rel.abs.bias$heights[min(posT$t)]  # height of current top strips

z.rel.abs.bias <- gtable_add_cols(z.rel.abs.bias, width, max(posR$r))  
z.rel.abs.bias <- gtable_add_rows(z.rel.abs.bias, height, min(posT$t)-1)

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, fill = "grey85")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 16, col = "grey10"))))

stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, fill = "grey85")),
  textGrob(labelT, gp = gpar(fontsize = 16, col = "grey10"))))

# Position the grobs in the gtable
z.rel.abs.bias <- gtable_add_grob(z.rel.abs.bias, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
z.rel.abs.bias <- gtable_add_grob(z.rel.abs.bias, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
z.rel.abs.bias <- gtable_add_cols(z.rel.abs.bias, unit(1/5, "line"), max(posR$r))
z.rel.abs.bias <- gtable_add_rows(z.rel.abs.bias, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(z.rel.abs.bias)

ggsave("plots/simulation_rel_abs_bias.png",plot = z.rel.abs.bias,scale=2)

# coverage
sim.results.coverage <- ggplot(data=results_long[results_long$variable=="CP",],
                           aes(x=factor(T0), y=value, colour=Estimator, group=forcats::fct_rev(Estimator)))  +   geom_line()  +
  facet_grid(shift_sc ~  NT, scales = "free", labeller=labeller)  +  xlab(TeX('Placebo $(a_i^{\\prime}/T)$')) + ylab("Coverage probability (%)") +
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

# Get the ggplot grob
z.coverage <- ggplotGrob(sim.results.coverage)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(z.coverage$layout, grepl("strip-r", name), select = t:r)
posT <- subset(z.coverage$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- z.coverage$widths[max(posR$r)]    # width of current right strips
height <- z.coverage$heights[min(posT$t)]  # height of current top strips

z.coverage <- gtable_add_cols(z.coverage, width, max(posR$r))  
z.coverage <- gtable_add_rows(z.coverage, height, min(posT$t)-1)

# Position the grobs in the gtable
z.coverage <- gtable_add_grob(z.coverage, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
z.coverage <- gtable_add_grob(z.coverage, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
z.coverage <- gtable_add_cols(z.coverage, unit(1/5, "line"), max(posR$r))
z.coverage <- gtable_add_rows(z.coverage, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(z.coverage)

ggsave("plots/simulation_coverage.png",plot = z.coverage,scale=2)

# boot_var

sim.results.boot.var <- ggplot(data=results_long[results_long$variable=="boot_var",],
                           aes(x=factor(T0), y=value, fill=Estimator))  + geom_boxplot(outlier.alpha = 0.3,outlier.size = 1, outlier.stroke = 0.1, lwd=0.25) +
  facet_grid(shift_sc ~  NT, scales = "free", labeller=labeller)  +  xlab(TeX('Placebo $(a_i^{\\prime}/T)$'))  + ylab("Bootstrap variance") +
  scale_fill_discrete(name = "Estimator:") +
  scale_y_continuous(breaks= pretty_breaks()) +
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

# Get the ggplot grob
z.boot.var <- ggplotGrob(sim.results.boot.var)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(z.boot.var$layout, grepl("strip-r", name), select = t:r)
posT <- subset(z.boot.var$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- z.boot.var$widths[max(posR$r)]    # width of current right strips
height <- z.boot.var$heights[min(posT$t)]  # height of current top strips

z.boot.var <- gtable_add_cols(z.boot.var, width, max(posR$r))  
z.boot.var <- gtable_add_rows(z.boot.var, height, min(posT$t)-1)

# Position the grobs in the gtable
z.boot.var <- gtable_add_grob(z.boot.var, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
z.boot.var <- gtable_add_grob(z.boot.var, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
z.boot.var <- gtable_add_cols(z.boot.var, unit(1/5, "line"), max(posR$r))
z.boot.var <- gtable_add_rows(z.boot.var, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(z.boot.var)

ggsave(paste0(output_dir,"simulation_boot_var.png",plot = z.boot.var,scale=2))