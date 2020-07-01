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

#covarflag <- c("","-covars")
covarflag <- c("")

for(o in outcome.vars){
  for(c in covarflag){
    print(c)
    print(o)
    
    ## Analysis 1: ST vs AT (retrospective, X=CBW) 
    
    print(paste0("Estimates for Analysis 1, outcome:",o,c))
    
    # get placebo results
    
    eastern.placebo.cbw <- readRDS(paste0("results/boot-trajectory-eastern-placebo-cbw",o,c,".rds"))
    
    swiss.placebo.cbw <- readRDS(paste0("results/boot-trajectory-eastern-placebo-cbw",o,c,".rds"))
    
    taus <- 1:length(eastern.placebo.cbw)
    
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
      labs(title="Retrospective prediction for later-treated", 
           subtitle=outcomes.labels[which(outcome.vars==o)],
           x="Randomization p-values",
           y=expression(delta)) + 
      geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
      scale_x_continuous(breaks=c(0,0.05,0.1,0.25,0.5,0.75,1), 
                         labels=c("0","0.05","0.10","0.25","0.50","0.75","1")) +
      scale_y_continuous(breaks=seq(3,15,3), 
                         labels=c("3",  "6", "9", "12","15"))+
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
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme_set(theme_bw() + theme(legend.key=element_blank(), legend.title=element_text(size=10))) + theme(plot.title = element_text(hjust = 0.5, size=14), plot.subtitle = element_text(hjust = 0.5, size=12)) + 
      theme(axis.text.y = element_text(size=8))
    
    ggsave(filename = paste0("plots/placebo-pvals-cbw-",o,c,".png"),plot = mc.placebo.plot)
 
    ## Analysis 2: ST vs NT (Prospective, X=LM) 
    
    print(paste0("Estimates for Analysis 2, outcome:",o,c))
    
    # get placebo results
    
    eastern.placebo.lm <- readRDS(paste0("results/boot-trajectory-eastern-placebo-lm",o,c,".rds"))
    
    swiss.placebo.lm <- readRDS(paste0("results/boot-trajectory-eastern-placebo-lm",o,c,".rds"))
    
    taus <- 1:length(eastern.placebo.lm)
    
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
      labs(title="Prospective prediction for later-treated", 
           subtitle=outcomes.labels[which(outcome.vars==o)],
           x="Randomization p-values",
           y= expression(delta)) + 
      geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
      scale_x_continuous(breaks=c(0,0.05,0.1,0.25,0.5,0.75,1), 
                         labels=c("0","0.05","0.10","0.25","0.50","0.75","1")) +
      scale_y_continuous(breaks=seq(3,15,3), 
                         labels=c("3",  "6", "9", "12","15"))+
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
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme_set(theme_bw() + theme(legend.key=element_blank(), legend.title=element_text(size=10))) + theme(plot.title = element_text(hjust = 0.5, size=14), plot.subtitle = element_text(hjust = 0.5, size=12)) + 
      theme(axis.text.y = element_text(size=8))
    
    ggsave(filename = paste0("plots/placebo-pvals-lm-",o,c,".png"),plot = mc.placebo.plot)
  }
}