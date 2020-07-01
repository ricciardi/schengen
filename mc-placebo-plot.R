# Plot p-values from placebo tests
library(ggplot2)
library(wesanderson)
library(reshape2)

outcome.vars <- c("N_CBWbord","CBWbord","CBWbordEMPL","empl","Thwusual","unempl","inact","seekdur_0","seekdur_1_2","seekdur_3more")
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
    print(o)
    
    ## Analysis 1: ST vs AT (retrospective, X=CBW) 
    
    print(paste0("Estimates for Analysis 1, outcome:",o,c))
    
    # get placebo results
    
    eastern.placebo.cbw <- readRDS(paste0("results/boot-trajectory-eastern-placebo-cbw-",o,c,".rds"))
    eastern.placebo.cbw.p <- sapply(1:length(eastern.placebo.cbw),  function(i) 1- ((1/length(eastern.placebo.cbw[[i]]$t) * sum(eastern.placebo.cbw[[i]]$t < eastern.placebo.cbw[[i]]$t0))))
    
    swiss.placebo.cbw <- readRDS(paste0("results/boot-trajectory-swiss-placebo-cbw-",o,c,".rds"))
    swiss.placebo.cbw.p <- sapply(1:length(swiss.placebo.cbw),  function(i) 1- ((1/length(swiss.placebo.cbw[[i]]$t) * sum(swiss.placebo.cbw[[i]]$t < swiss.placebo.cbw[[i]]$t0))))
    
    taus <- 1:length(eastern.placebo.cbw)
    
    p.values <- data.frame("Eastern"=c(eastern.placebo.cbw.p),
                           "Swiss"=c(swiss.placebo.cbw.p),
                           "tau"=taus)
    
    p.values.m <-melt(p.values,id.vars = c("tau"))
    
    # Plot
    mc.placebo.plot <- ggplot(p.values.m, aes(x=value, y=tau)) + 
      geom_point(stat='identity', aes(col=variable,shape=factor(variable)), size=3, alpha=0.5)  +
      labs(title="Placebo test: retrospective prediction for later-treated", 
           subtitle=outcomes.labels[which(outcome.vars==o)],
           x="Block bootstrap p-values",
           y=expression(rho)) + 
      geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
      scale_x_continuous(breaks=c(0,0.05,0.1,0.25,0.5,0.75,1), 
                         labels=c("0","0.05","0.10","0.25","0.50","0.75","1")) +
      coord_flip() +
      scale_shape_manual(name="Cluster", values = c("Eastern" = 2,
                                                    "Swiss" = 4),
                         labels=c("Eastern", "Swiss")) +
      scale_colour_manual(name="Cluster", values = c(  "Eastern" = wes_palette("Darjeeling1")[1],
                                                       "Swiss" = wes_palette("Darjeeling1")[2]),
                          labels=c("Eastern", "Swiss")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme_set(theme_bw() + theme(legend.key=element_blank(), legend.title=element_text(size=10))) + theme(plot.title = element_text(hjust = 0.5, size=14), plot.subtitle = element_text(hjust = 0.5, size=12)) + 
      theme(axis.text.y = element_text(size=8))
    
    ggsave(filename = paste0("plots/placebo-pvals-cbw-",o,c,".png"),plot = mc.placebo.plot)
    
    ## Analysis 2: ST vs NT (Prospective, X=LM) 
    
    if(o %in% c("N_lmbord","lmbord","lmbordEMPL")) next
    
    print(paste0("Estimates for Analysis 2, outcome:",o,c))
    
    # get placebo results
    
    eastern.placebo.lm <- readRDS(paste0("results/boot-trajectory-eastern-placebo-lm-",o,c,".rds"))
    eastern.placebo.lm.p <- sapply(1:length(eastern.placebo.lm),  function(i) 1- ((1/length(eastern.placebo.lm[[i]]$t) * sum(eastern.placebo.lm[[i]]$t < eastern.placebo.lm[[i]]$t0))))
    
    swiss.placebo.lm <- readRDS(paste0("results/boot-trajectory-swiss-placebo-lm-",o,c,".rds"))
    swiss.placebo.lm.p <- sapply(1:length(swiss.placebo.lm),  function(i) 1- ((1/length(swiss.placebo.lm[[i]]$t) * sum(swiss.placebo.lm[[i]]$t < swiss.placebo.lm[[i]]$t0))))
    
    taus <- 1:length(eastern.placebo.lm)
    
    p.values <- data.frame("Eastern"=c(eastern.placebo.lm.p),
                           "Swiss"=c(swiss.placebo.lm.p),
                           "tau"=taus)
    
    p.values.m <-melt(p.values,id.vars = c("tau"))
    
    # Plot
    mc.placebo.plot <- ggplot(p.values.m, aes(x=value, y=tau)) + 
      geom_point(stat='identity', aes(col=variable,shape=factor(variable)), size=3, alpha=0.5)  +
      labs(title="Placebo test: prospective prediction for later-treated", 
           subtitle=outcomes.labels[which(outcome.vars==o)],
           x="Block bootstrap p-values",
           y=expression(rho)) + 
      geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
      scale_x_continuous(breaks=c(0,0.05,0.1,0.25,0.5,0.75,1), 
                         labels=c("0","0.05","0.10","0.25","0.50","0.75","1")) +
      coord_flip() +
      scale_shape_manual(name="Cluster", values = c("Eastern" = 2,
                                                    "Swiss" = 4),
                         labels=c("Eastern", "Swiss")) +
      scale_colour_manual(name="Cluster", values = c(  "Eastern" = wes_palette("Darjeeling1")[1],
                                                       "Swiss" = wes_palette("Darjeeling1")[2]),
                          labels=c("Eastern", "Swiss")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme_set(theme_bw() + theme(legend.key=element_blank(), legend.title=element_text(size=10))) + theme(plot.title = element_text(hjust = 0.5, size=14), plot.subtitle = element_text(hjust = 0.5, size=12)) + 
      theme(axis.text.y = element_text(size=8))
    
    ggsave(filename = paste0("plots/placebo-pvals-lm-",o,c,".png"),plot = mc.placebo.plot)
  }
}