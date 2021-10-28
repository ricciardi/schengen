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

# Load results data

filenames <- c(list.files(path="results/simulation", pattern = ".rds", full.names = TRUE))

results <- list() # structure is: [[filename]][[metric]]
for(f in filenames){
  print(f)
  result.matrix <- readRDS(f)
  n <- nrow(result.matrix )
  bias <- matrix(NA, n, 3)
  CP <- matrix(NA, n, 3)
  RMSE <- matrix(NA, n, 3)
  for(i in c('est_mc_plain_bias','est_mc_weights_bias','est_mc_weights_covars_bias')){
    bias[,which(c('est_mc_plain_bias','est_mc_weights_bias','est_mc_weights_covars_bias')==i)] <- unlist(result.matrix[,i])
  }
  for(i in c('est_mc_plain_cp','est_mc_weights_cp','est_mc_weights_covars_cp')){
    CP[,which(c('est_mc_plain_cp','est_mc_weights_cp','est_mc_weights_covars_cp')==i)] <- unlist(result.matrix[,i])
  }
  for(i in c('est_mc_plain_RMSE','est_mc_weights_RMSE','est_mc_weights_covars_RMSE')){
    RMSE[,which(c('est_mc_plain_RMSE','est_mc_weights_RMSE','est_mc_weights_covars_RMSE')==i)] <- unlist(result.matrix[,i])
  }
  fr_obs <- unlist(result.matrix[,"fr_obs"])
  results[[f]] <- list("bias"=bias,"CP"=CP,"RMSE"=RMSE,"fr_obs"=fr_obs,"n"=n)
}

# Create New lists
# structure is: [[estimator]][[filename]]

bias <- list()
for(i in 1:length(c('est_mc_plain_bias','est_mc_weights_bias','est_mc_weights_covars_bias'))){
  bias[[i]] <- lapply(1:length(filenames), function(f) results[[f]]$bias[,i])
}

CP <- list()
for(i in 1:length(c('est_mc_plain_cp','est_mc_weights_cp','est_mc_weights_covars_cp'))){
  CP[[i]] <- lapply(1:length(filenames), function(f) results[[f]]$CP[,i])
}

RMSE <- list()
for(i in 1:length(c('est_mc_plain_RMSE','est_mc_weights_RMSE','est_mc_weights_covars_RMSE'))){
  RMSE[[i]] <- lapply(1:length(filenames), function(f) results[[f]]$RMSE[,i])
}

RMSE <- list()
for(i in 1:length(c('est_mc_plain_RMSE','est_mc_weights_RMSE','est_mc_weights_covars_RMSE'))){
  RMSE[[i]] <- lapply(1:length(filenames), function(f) results[[f]]$RMSE[,i])
}

#fr_obs <- sapply(1:length(filenames), function(f) results[[f]]$fr_obs)


# Create dataframe for plot
results.df <- data.frame("abs.bias"=abs(unlist(bias)),
                         "Coverage"=unlist(CP),
                         "RMSE"=unlist(RMSE),
                         "fr_obs"= rep(sapply(1:length(filenames), function(f) results[[f]]$fr_obs), 3),
                         "Estimator"=c(rep("MC (Athey et al.)",length.out=length(c(unlist(CP[[1]])))), rep("MC (weighted loss)",length.out=length(c(unlist(CP[[2]])))),
                         rep("MC (weighted loss + covariates)",length.out=length(c(unlist(CP[[3]]))))),
                         "filename"=c(rep(unlist(sapply(1:length(filenames), function(i) rep(filenames[i], length.out=n))), 3)))

results.df$NT <- NA
results.df$R <- NA
results.df$noise_sc <- NA
results.df$effect_size <- NA
  
for(s in c("N_40_T_40","N_50_T_50","N_60_T_60")){
  if(length(grep(s, results.df$filename))>0){
    print(s)
    if(s=="N_40_T_40"){
      NT <- 40**2
    }
    if(s=="N_50_T_50"){
      NT <- 50**2
    }
    if(s=="N_60_T_60"){
      NT <- 60**2
    }
    results.df[grep(s, results.df$filename),]$NT <- as.numeric(NT)
  }
}

for(s in c("R_20","R_30","R_40")){
  if(length(grep(s, results.df$filename))>0){
    print(s)
    if(s=="R_20"){
      R <- 20
    }
    if(s=="R_30"){
      R <- 30
    }
    if(s=="R_40"){
      R <- 40
    }
    results.df[grep(s, results.df$filename),]$R <- R
  }
}

for(s in c("noise_sc_0.1","noise_sc_0.2","noise_sc_0.3")){
  if(length(grep(s, results.df$filename))>0){
    print(s)
    if(s=="noise_sc_0.1"){
      noise_sc <- 0.1
    }
    if(s=="noise_sc_0.2"){
      noise_sc <- 0.2
    }
    if(s=="noise_sc_0.3"){
      noise_sc <- 0.3
    }
    results.df[grep(s, results.df$filename),]$noise_sc  <- noise_sc
  }
}

for(s in c("effect_size_0.001","effect_size_0.005","effect_size_0.01")){
  if(length(grep(s, results.df$filename))>0){
    print(s)
    if(s=="effect_size_0.001"){
      effect_size <- 0.001
    }
    if(s=="effect_size_0.005"){
      effect_size <- 0.005
    }
    if(s=="effect_size_0.01"){
      effect_size <- 0.01
    }
    results.df[grep(s, results.df$filename),]$effect_size  <- effect_size
  }
}

# create coverage rate variable

results.df <- results.df %>% 
  group_by(Estimator,NT,R,noise_sc,effect_size) %>% 
  mutate(CP = mean(Coverage)) 

# reshape and plot
results.df$id <- with(results.df, paste(NT,R,noise_sc,effect_size, sep = "_"))
results_long <- reshape2::melt(results.df[!colnames(results.df) %in% c("id","filename")], id.vars=c("Estimator","NT","R","noise_sc","effect_size","fr_obs"))  # convert to long format

# bias (NxT)
sim.results.bias <- ggplot(data=results_long[results_long$variable=="abs.bias" & results_long$NT==1600,],
                           aes(x=factor(R), y=value, fill=forcats::fct_rev(Estimator)))  + geom_boxplot(outlier.alpha = 0.6,outlier.size = 1.2, outlier.stroke = 0.2, lwd=0.25) +
  facet_grid(factor(noise_sc) ~  factor(effect_size), scales = "fixed")  +  xlab("Rank (R)") + ylab("Absolute bias") + ggtitle(TeX("Absolute bias, $N \\times T = 1600$")) +
  scale_fill_discrete(name = "Estimator:") +
  theme(legend.position="bottom") +   theme(plot.title = element_text(hjust = 0.5, family="serif", size=16)) +
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
z.bias <- ggplotGrob(sim.results.bias)

# Labels 
labelR <- "Degree of confounding"
labelT <- "Effect size"

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(z.bias$layout, grepl("strip-r", name), select = t:r)
posT <- subset(z.bias$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- z.bias$widths[max(posR$r)]    # width of current right strips
height <- z.bias$heights[min(posT$t)]  # height of current top strips

z.bias <- gtable_add_cols(z.bias, width, max(posR$r))  
z.bias <- gtable_add_rows(z.bias, height, min(posT$t)-1)

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, fill = "grey85")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 16, col = "grey10"))))

stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, fill = "grey85")),
  textGrob(labelT, gp = gpar(fontsize = 16, col = "grey10"))))

# Position the grobs in the gtable
z.bias <- gtable_add_grob(z.bias, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
z.bias <- gtable_add_grob(z.bias, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
z.bias <- gtable_add_cols(z.bias, unit(1/5, "line"), max(posR$r))
z.bias <- gtable_add_rows(z.bias, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(z.bias)

ggsave("plots/simulation_bias.png",plot = z.bias,scale=2)

# coverage
sim.results.coverage <- ggplot(data=results_long[results_long$variable=="CP" & results_long$NT==1600,],
                           aes(x=factor(R), y=value, colour=forcats::fct_rev(Estimator), group=forcats::fct_rev(Estimator)))  +   geom_line()  +
  facet_grid(factor(noise_sc) ~  factor(effect_size), scales = "fixed")  +  xlab("Rank (R)") + ylab("Coverage (%)") + ggtitle(TeX("Coverage, $N \\times T = 1600$")) + 
  scale_colour_discrete(name = "Estimator:") +
  geom_hline(yintercept = 0.95, linetype="dotted")+
  theme(legend.position="bottom") +   theme(plot.title = element_text(hjust = 0.5, family="serif", size=16)) +
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

# RMSE
  
sim.results.RMSE <- ggplot(data=results_long[results_long$variable=="RMSE" & results_long$NT==1600,],
                           aes(x=factor(R), y=value, fill=forcats::fct_rev(Estimator)))  + geom_boxplot(outlier.alpha = 0.6,outlier.size = 1.2, outlier.stroke = 0.2, lwd=0.25) +
  facet_grid(factor(noise_sc) ~  factor(effect_size), scales = "fixed")  +  xlab("Rank (R)")  + ylab("RMSE") + ggtitle(TeX("RMSE, $N \\times T = 1600$")) +
  scale_fill_discrete(name = "Estimator:") +
  theme(legend.position="bottom") +   theme(plot.title = element_text(hjust = 0.5, family="serif", size=16)) +
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
z.RMSE <- ggplotGrob(sim.results.RMSE)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(z.RMSE$layout, grepl("strip-r", name), select = t:r)
posT <- subset(z.RMSE$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- z.RMSE$widths[max(posR$r)]    # width of current right strips
height <- z.RMSE$heights[min(posT$t)]  # height of current top strips

z.RMSE <- gtable_add_cols(z.RMSE, width, max(posR$r))  
z.RMSE <- gtable_add_rows(z.RMSE, height, min(posT$t)-1)

# Position the grobs in the gtable
z.RMSE <- gtable_add_grob(z.RMSE, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
z.RMSE <- gtable_add_grob(z.RMSE, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
z.RMSE <- gtable_add_cols(z.RMSE, unit(1/5, "line"), max(posR$r))
z.RMSE <- gtable_add_rows(z.RMSE, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(z.RMSE)

ggsave("plots/simulation_RMSE.png",plot = z.RMSE,scale=2)