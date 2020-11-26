TsPlot <- function(df, main = "",y.title,vline,vline2,breaks,labels,hline) {
  library(ggplot2)
  library(zoo)
  library(scales)
  library(wesanderson)
  
  gg.xts <- ggplot(df, aes(x = quarter)) +
    
    # panel layout
    facet_grid(series~., scales = "free_y", space = "fixed", shrink = TRUE, drop = TRUE, labeller=label_parsed) + 
    
    theme(strip.text= element_text(size = 16, family = "serif", face='bold'), strip.background = element_blank()) +
    
    # line colours
    geom_line(data = subset(df, variable == "observed.eastern"), aes(y = value, colour = "observed.eastern", linetype="observed.eastern"), show.legend = TRUE, size=1) +
    geom_point(data = subset(df, variable == "predicted.eastern"), aes(y = value, colour = "predicted.eastern", shape="predicted.eastern"), show.legend = FALSE, alpha=0.8, size=2.5) +
    geom_line(data = subset(df, variable == "pointwise.eastern"), aes(y = value, colour = "predicted.eastern", linetype="predicted.eastern"), show.legend = FALSE, size=1, na.rm=TRUE) +
    
    geom_line(data = subset(df, variable == "observed.swiss"), aes(y = value, colour = "observed.swiss", linetype="observed.swiss"), show.legend = TRUE, size=1) +
    geom_point(data = subset(df, variable == "predicted.swiss"), aes(y = value, colour = "predicted.swiss", shape="predicted.swiss"), show.legend = FALSE, alpha=0.8, size=2.5) +
    geom_line(data = subset(df, variable == "pointwise.swiss"), aes(y = value, colour = "predicted.swiss", linetype="predicted.swiss"), show.legend = FALSE, size=1, na.rm = TRUE) +
    
    geom_line(data = subset(df, variable == "observed.control"), aes(y = value, colour = "observed.control", linetype="observed.control"), show.legend = TRUE, size=1) +
    
    # intervals
    geom_ribbon(data = subset(df, variable == "pointwise.eastern"), aes(ymin = lower, ymax=upper, colour="predicted.eastern"), alpha=.1, size=0.2, show.legend = FALSE) +
    geom_ribbon(data = subset(df, variable == "pointwise.swiss"), aes(ymin = lower, ymax=upper, colour="predicted.swiss"), alpha=.1, size=0.2, show.legend = FALSE) +
    
    # horizontal line to indicate zero values
    geom_hline(aes(yintercept = hline), size = 0.5, colour = "black") +
    
    # main y-axis title
    ylab(y.title) +
    
    # x-axis title
    xlab("YearQuarter") +
    
    # main chart title
    ggtitle(main)
  
  # vertical line to indicate intervention
  
  intervention <- geom_vline(xintercept=vline, linetype=3)
  intervention2 <- geom_vline(xintercept=vline2, linetype=2)
  
  # horizontal ticks
  
  ticks <- scale_x_continuous(breaks=breaks,
                            labels=labels)
  
  # annotation text
  
  # ann_text <- data.frame(year = c(20074, 20084), value=0.001,
  #                        series = factor("Time-series", levels = c("Time-series", "Per-period effect")),
  #                        lab = c("Pre","Post"))
  
  # legend 
  lines <- scale_linetype_manual(name="", values = c("observed.control" = "dashed",
                                                     "observed.eastern" = "solid",
                                                     "observed.swiss" = "solid",
                                                     "predicted.eastern" = "dotted",
                                                     "predicted.swiss" = "dotted"),
                                 labels=c("Observed Always-treated", "Observed Eastern", "Observed Swiss",
                                          "Predicted Eastern", "Predicted Swiss")) 
  
  shapes <- scale_shape_manual(name="", values = c("predicted.eastern" = 1,
                                                   "predicted.swiss" = 2),
                               labels=c("Predicted Eastern", "Predicted Swiss")) 
  
  colours <-     scale_colour_manual(name="", values = c(  "observed.control" = wes_palette("Darjeeling1")[1],
                                                           "observed.eastern" = wes_palette("Darjeeling1")[5], 
                                                           "observed.swiss" = wes_palette("Darjeeling1")[4],
                                                           "predicted.eastern" = wes_palette("Darjeeling1")[5],
                                                           "predicted.swiss" = wes_palette("Darjeeling1")[4]),
                                     labels=c("Observed Always-treated", "Observed Eastern", "Observed Swiss",
                                              "Predicted Eastern", "Predicted Swiss")) 
  gg.xts <- gg.xts +
    intervention + intervention2 +
    ticks + 
    theme( legend.title = element_blank()
           , plot.title = element_text(hjust = 0.5, size = 16)
           , legend.justification = c(0.98,0.25)
           , legend.position = c(0.99,0.25)
           #  , legend.position = "top"
           , legend.background = element_rect(fill="transparent")
           , axis.text=element_text(size=10)
           , axis.title.x=element_text(size = 14, margin = margin(t = 20, r = 0, b = 0, l = 0))
           , axis.title.y=element_text(size = 14, margin = margin(t = 0, r = 20, b = 0, l = 0))
           , legend.text=element_text(size=14, family = "serif")
           , legend.box = "vertical"
           , legend.key = element_blank()
    ) +
    #+ geom_text(data = ann_text,aes(y = value, label =lab), family="serif", fontface="italic",  size=6) +
    #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #          panel.background = element_blank(), axis.line = element_line(colour = "black")) + # rm background
    colours+ 
    lines +
    shapes + 
    theme(legend.key.width=unit(4,"line"))
  
  return(gg.xts)
}