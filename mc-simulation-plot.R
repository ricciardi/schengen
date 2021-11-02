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

n.estimators <- 5

# Load results data

filenames <- c(list.files(path="mc-simulations", pattern = ".rds", full.names = TRUE))
filenames <- filenames[grep("2000",filenames)]

bias.vars <- c('est_mc_plain_bias','est_mc_weights_bias','est_mc_weights_covars_bias','est_model_ADH_bias','est_model_DID_bias')
cp.vars <- c('est_mc_plain_cp','est_mc_weights_cp','est_mc_weights_covars_cp','est_model_ADH_cp','est_model_DID_cp')
RMSE.vars <- c('est_mc_plain_RMSE','est_mc_weights_RMSE','est_mc_weights_covars_RMSE','est_model_ADH_RMSE','est_model_DID_RMSE')
ciw.vars <- c('est_mc_plain_ciw','est_mc_weights_ciw','est_mc_weights_covars_ciw','est_model_ADH_ciw','est_model_DID_ciw')

results <- list() # structure is: [[filename]][[metric]]
for(f in filenames){
  print(f)
  result.matrix <- readRDS(f)
  n <- nrow(result.matrix )
  bias <- matrix(NA, n, n.estimators)
  CP <- matrix(NA, n, n.estimators)
  RMSE <- matrix(NA, n, n.estimators)
  CIW <- matrix(NA, n, n.estimators)
  for(i in bias.vars){
    bias[,which(bias.vars==i)] <- unlist(result.matrix[,i])
  }
  for(i in cp.vars){
    CP[,which(cp.vars==i)] <- unlist(result.matrix[,i])
  }
  for(i in RMSE.vars){
    RMSE[,which(RMSE.vars==i)] <- unlist(result.matrix[,i])
  }
  for(i in ciw.vars){
    CIW[,which(ciw.vars==i)] <- unlist(result.matrix[,i])
  }
  fr_obs <- unlist(result.matrix[,"fr_obs"])
  results[[f]] <- list("bias"=bias,"CP"=CP,"RMSE"=RMSE,"CIW"=CIW,"fr_obs"=fr_obs,"n"=n)
}

# Create New lists
# structure is: [[estimator]][[filename]]

bias <- list()
for(i in 1:length(bias.vars)){
  bias[[i]] <- lapply(1:length(filenames), function(f) results[[f]]$bias[,i])
}

CP <- list()
for(i in 1:length(cp.vars)){
  CP[[i]] <- lapply(1:length(filenames), function(f) results[[f]]$CP[,i])
}

RMSE <- list()
for(i in 1:length(RMSE.vars)){
  RMSE[[i]] <- lapply(1:length(filenames), function(f) results[[f]]$RMSE[,i])
}

CIW <- list()
for(i in 1:length(ciw.vars)){
  CIW[[i]] <- lapply(1:length(filenames), function(f) results[[f]]$CIW[,i])
}

# Create dataframe for plot
results.df <- data.frame("abs.bias"=abs(as.numeric(unlist(bias))),
                         "Coverage"=as.numeric(unlist(CP)),
                         "RMSE"=as.numeric(unlist(RMSE)),
                         "CIW"=as.numeric(unlist(CIW)),
                         "fr_obs"= rep(sapply(1:length(filenames), function(f) results[[f]]$fr_obs), n.estimators),
                         "Estimator"=c(rep("MC (Athey et al.)",length.out=length(c(unlist(CP[[1]])))), 
                                       rep("MC (weighted)",length.out=length(c(unlist(CP[[2]])))),
                         rep("MC (weighted + covars.)",length.out=length(c(unlist(CP[[3]])))),
                         rep("SCM",length.out=length(c(unlist(CP[[4]])))),
                         rep("DID",length.out=length(c(unlist(CP[[5]]))))),
                         "filename"=c(rep(unlist(sapply(1:length(filenames), function(i) rep(filenames[i], length.out=n))), n.estimators)))

results.df$NT <- NA
results.df$R <- NA
results.df$noise_sc <- NA
  
for(s in c("N_40_T_40","N_60_T_60","N_80_T_80")){
  if(length(grep(s, results.df$filename))>0){
    print(s)
    if(s=="N_40_T_40"){
      NT <- "1600"
    }
    if(s=="N_60_T_60"){
      NT <- "3600"
    }
    if(s=="N_80_T_80"){
      NT <- "6400"
    }
    results.df[grep(s, results.df$filename),]$NT <- NT
  }
}

for(s in c("R_10","R_20","R_40")){
  if(length(grep(s, results.df$filename))>0){
    print(s)
    if(s=="R_10"){
      R <- "10"
    }
    if(s=="R_20"){
      R <- "20"
    }
    if(s=="R_40"){
      R <- "40"
    }
    results.df[grep(s, results.df$filename),]$R <- R
  }
}

for(s in c("noise_sc_0.01","noise_sc_0.1","noise_sc_0.2")){
  if(length(grep(s, results.df$filename))>0){
    print(s)
    if(s=="noise_sc_0.01"){
      noise_sc <- "0.01"
    }
    if(s=="noise_sc_0.1"){
      noise_sc <- "0.1"
    }
    if(s=="noise_sc_0.2"){
      noise_sc <- "0.2"
    }
    results.df[grep(s, results.df$filename),]$noise_sc  <- noise_sc
  }
}

# create coverage rate variable

results.df <- results.df %>% 
  group_by(Estimator,NT,R,noise_sc) %>% 
  mutate(CP = mean(Coverage)) 

# reshape and plot
results.df$id <- with(results.df, paste(NT,R,noise_sc, sep = "_"))
results_long <- reshape2::melt(results.df[!colnames(results.df) %in% c("id","filename")], id.vars=c("Estimator","NT","R","noise_sc","fr_obs"))  # convert to long format

variable_names <- list(
  '1600'= TeX("$40 \\times 40"),
  '3600'= TeX("$60 \\times 60"),
  '6400'= TeX("$60 \\times 60"),
  '0.01'= '0.01',
  '0.1'= '0.1',
  '0.2'= '0.2') 

labeller <- function(variable,value){
  return(variable_names[value])
}
# bias (NxT)
sim.results.bias <- ggplot(data=results_long[results_long$variable=="abs.bias",],
                           aes(x=factor(R), y=value, fill=Estimator))  + geom_boxplot(outlier.alpha = 0.3,outlier.size = 1, outlier.stroke = 0.1, lwd=0.25) +
  facet_grid(noise_sc ~  NT, scales = "free", labeller=labeller)  +  xlab("Rank (R)") + ylab("Absolute bias") +
  scale_fill_discrete(name = "Estimator:") +
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
z.bias <- ggplotGrob(sim.results.bias)

# Labels 
labelR <- "Noise scale"
labelT <- TeX("$N \\times T")

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
sim.results.coverage <- ggplot(data=results_long[results_long$variable=="CP",],
                           aes(x=factor(R), y=value, colour=Estimator, group=forcats::fct_rev(Estimator)))  +   geom_line()  +
  facet_grid(noise_sc ~  NT, scales = "free", labeller=labeller)  +  xlab("Rank (R)") + ylab("Coverage probability (%)") +
  scale_colour_discrete(name = "Estimator:") +
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

# RMSE
  
sim.results.RMSE <- ggplot(data=results_long[results_long$variable=="RMSE",],
                           aes(x=factor(R), y=value, fill=forcats::fct_rev(Estimator)))  + geom_boxplot(outlier.alpha = 0.3,outlier.size = 1, outlier.stroke = 0.1, lwd=0.25) +
  facet_grid(noise_sc ~  NT, scales = "free", labeller=labeller)  +  xlab("Rank (R)")  + ylab("RMSE") +
  scale_fill_discrete(name = "Estimator:") +
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

# CIW

sim.results.CIW <- ggplot(data=results_long[results_long$variable=="CIW",],
                           aes(x=factor(R), y=value, fill=forcats::fct_rev(Estimator)))  + geom_boxplot(outlier.alpha = 0.3,outlier.size = 1, outlier.stroke = 0.1, lwd=0.25) +
  facet_grid(noise_sc ~  NT, scales = "free", labeller=labeller)  +  xlab("Rank (R)")  + ylab("Confidence interval width") +
  scale_fill_discrete(name = "Estimator:") +
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
z.CIW <- ggplotGrob(sim.results.CIW)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(z.CIW$layout, grepl("strip-r", name), select = t:r)
posT <- subset(z.CIW$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- z.CIW$widths[max(posR$r)]    # width of current right strips
height <- z.CIW$heights[min(posT$t)]  # height of current top strips

z.CIW <- gtable_add_cols(z.CIW, width, max(posR$r))  
z.CIW <- gtable_add_rows(z.CIW, height, min(posT$t)-1)

# Position the grobs in the gtable
z.CIW <- gtable_add_grob(z.CIW, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
z.CIW <- gtable_add_grob(z.CIW, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
z.CIW <- gtable_add_cols(z.CIW, unit(1/5, "line"), max(posR$r))
z.CIW <- gtable_add_rows(z.CIW, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(z.CIW)

ggsave("plots/simulation_CIW.png",plot = z.CIW,scale=2)