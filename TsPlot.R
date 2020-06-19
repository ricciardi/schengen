TsPlot <- function(df, main = "",y.title,vline,breaks,labels) {
  library(ggplot2)
  library(zoo)
  library(scales)
  library(wesanderson)
  
  gg.xts <- ggplot(df, aes(x = year)) +
    
    # panel layout
    facet_grid(series~., scales = "free_y", space = "fixed", shrink = TRUE, drop = TRUE, labeller=label_parsed) + 
    
    theme(strip.text= element_text(size = 16, family = "serif", face='bold'), strip.background = element_blank()) +
    
    # line colours
    geom_line(data = subset(df, variable == "observed.treated"), aes(y = value, colour = "observed.treated", linetype="observed.treated"), show.legend = TRUE, size=1) +
    geom_line(data = subset(df, variable == "predicted.treated"), aes(y = value, colour = "predicted.treated", linetype="predicted.treated"), show.legend = FALSE, size=1) +
    geom_line(data = subset(df, variable == "pointwise.treated"), aes(y = value, colour = "predicted.treated", linetype="predicted.treated"), show.legend = FALSE, size=1) +

    geom_line(data = subset(df, variable == "observed.control"), aes(y = value, colour = "observed.control", linetype="observed.control"), show.legend = TRUE, size=1) +
    
    # intervals
    geom_ribbon(data = subset(df, variable == "pointwise.treated"), aes(ymin = lower, ymax=upper, colour="predicted.treated"), alpha=.1, size=0.5, show.legend = FALSE) +

    # horizontal line to indicate zero values
    geom_hline(yintercept = 0, size = 0.5, colour = "black") +
    
    # main y-axis title
    ylab(y.title) +
    
    # x-axis title
    xlab("\nTime") +
  
    # main chart title
    ggtitle(main)
  
  # vertical line to indicate intervention
  
  intervention <- geom_vline(xintercept=vline, linetype=2)
  
  # horizontal ticks
  
  ticks <- scale_x_continuous(breaks=breaks,
                           labels=labels)
  
  # annotation text
  
  # ann_text <- data.frame(year = c(20074, 20084), value=0.001,
  #                        series = factor("Time-series", levels = c("Time-series", "Per-period effect")),
  #                        lab = c("pre-period","post-period"))
  
  if(rev){
  # legend 
    lines <- scale_linetype_manual(name="", values = c("observed.control" = "dashed",
                                                       "observed.treated" = "solid",
                                                       "predicted.treated" = "dotted"),
                                   labels=c("Observed always-treated", "Observed later-treated", 
                                            "Predicted later-treated")) 
    colours <-     scale_colour_manual(name="", values = c(  "observed.control" = wes_palette("Darjeeling1")[1],
                                                             "observed.treated" = wes_palette("Darjeeling1")[5], 
                                                             "predicted.treated" = wes_palette("Darjeeling1")[5]),
                                       labels=c("Observed always-treated", "Observed later-treated", 
                                                "Predicted later-treated")) 
  }else{
    lines <- scale_linetype_manual(name="", values = c("observed.control" = "dashed",
                                                       "observed.treated" = "solid",
                                                       "predicted.treated" = "dotted"),
                                   labels=c("Observed never-treated", "Observed later-treated", 
                                            "Predicted later-treated")) 
    colours <-     scale_colour_manual(name="", values = c(  "observed.control" = wes_palette("Darjeeling1")[1],
                                                             "observed.treated" = wes_palette("Darjeeling1")[5], 
                                                             "predicted.treated" = wes_palette("Darjeeling1")[5]),
                                       labels=c("Observed never-treated", "Observed later-treated", 
                                                "Predicted later-treated")) 
  }
  gg.xts <- gg.xts +
    intervention +
    ticks + 
    theme( legend.title = element_blank()
           , plot.title = element_text(hjust = 0.5)
           , legend.justification = c(1,0)
            , legend.position = "bottom"
           , legend.background = element_rect()
           , axis.text=element_text(size=14)
           , axis.title.x=element_text(size = 16)
           , axis.title.y=element_text(size = 16)
           , legend.text=element_text(size=14, family = "serif")
           , legend.box = "horizontal" 
    ) +
   #+ geom_text(data = ann_text,aes(y = value, label =lab), family="serif", fontface="italic",  size=6) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + # rm background
    colours+ 
    lines +
    theme(legend.key.width=unit(4,"line")) 
  return(gg.xts)
}