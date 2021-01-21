library(ggplot2)
library(latex2exp)
library(dplyr)
library(scales)

outcome.vars <- c("CBWbord","CBWbordEMPL")
outcomes.labels <- c("% working in border region",
                     "% working in border region,\n conditional on employment")

sim.labels <- c("Staggered adoption","Simultaneous adoption")

for(o in outcome.vars){
  for(i in c(0,1)){
    print(o)
    print(i)
    
    if(o=="CBWbord"){
      load(paste0("results/CBWbord_N_53_T_36_numruns_100_num_treated_27_simultaneuous_",i,"_covars_FALSE.rds"))
    }
    if(o=="CBWbordEMPL"){
      load(paste0("results/CBWbordEMPL_N_53_T_36_numruns_100_num_treated_27_simultaneuous_",i,"_covars_FALSE.rds"))
    }
    
    df1 <- df1 %>% group_by(x) %>% mutate(y = y,
                                          lb= lb,
                                          ub = ub)
    df1$x <- round(df1$x,2)
    levels(df1$Method) <- c("DID", "MC-NNM", "SCM")
    
    schengen <- ggplot(data = df1, aes(x, y, color = Method, shape = Method)) +
      geom_point(size = 5, position=position_dodge(width=0.1)) +
      geom_errorbar(
        aes(ymin = lb, ymax = ub),
        width = 0.1,
        linetype = "solid",
        position=position_dodge(width=0.1)) +
      scale_shape_manual(values=c(1:3)) +
      scale_x_continuous(breaks=c(unique(df1$x)[1],unique(df1$x)[2],unique(df1$x)[3],unique(df1$x)[4],unique(df1$x)[5])) +
      scale_y_continuous(breaks= pretty_breaks()) +
      theme_bw() +
      xlab(TeX('Placebo $(T_0/T)$')) +
      ylab("Average RMSE") +
      theme(axis.title=element_text(family="serif", size=16)) +
      theme(axis.text=element_text(family="serif", size=14)) +
      theme(legend.text=element_text(family="serif", size = 12)) +
      theme(legend.title=element_text(family="serif", size = 12)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l =0))) +
      theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l =0))) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))  # rm background
    
    ggsave(paste0("plots/",o,"-sim-",i,".png"), schengen + theme(legend.position = "none"), scale=1.25)
    ggsave(paste0("plots/",o,"-sim-",i,"slides.png"), schengen + ggtitle(paste0("Placebo test in retrospective setting"),subtitle=outcomes.labels[which(outcome.vars==o)]) + theme(plot.title = element_text(family="serif", size=16, hjust = 0.5), plot.subtitle = element_text(family="serif", size=12, hjust = 0.5)), scale=1.25)
  }
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_color_hue(3)