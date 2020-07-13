# Plot CIs from placebo tests
library(ggplot2)
library(wesanderson)
library(reshape2)

outcome.vars <- c("CBWbordEMPL","empl","Thwusual","unempl","seekdur_3more")
outcomes.labels <- c("% working in border region", # order of covar name
                     "Employment rate",
                     "% unemployed for\n < 1 year",
                     "Average total\n working hours",
                     "Unemployment rate")

cf<- ("")
outcome.pvals <- lapply(outcome.vars, function(o){lapply
  
  #Analysis 1: ST vs AT (retrospective, X=CBW) 
  
  print(paste0("Estimates for Analysis 1, outcome:",o,cf))
  
  # placebo data
  outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  # Use pre-treatment
  outcomes.cbw.placebo <- outcomes.cbw
  outcomes.cbw.placebo$mask <- outcomes.cbw$mask[,1:(which(colnames(outcomes.cbw$mask)=="20072")-1)]
  outcomes.cbw.placebo$mask[outcomes.cbw.placebo$mask>0] <- 0
  outcomes.cbw.placebo$mask[rownames(outcomes.cbw.placebo$mask)%in%outcomes.cbw.placebo$eastern,][,1:(which(colnames(outcomes.cbw.placebo$mask)=="20052"))] <- 1
  outcomes.cbw.placebo$mask[rownames(outcomes.cbw.placebo$mask)%in%outcomes.cbw.placebo$swiss,][,1:(which(colnames(outcomes.cbw.placebo$mask)=="20054"))] <- 1
  outcomes.cbw.placebo$M <- outcomes.cbw$M[,1:(which(colnames(outcomes.cbw$mask)=="20072")-1)]
  outcomes.cbw.placebo$W <- outcomes.cbw$W[,1:(which(colnames(outcomes.cbw$mask)=="20072")-1)]
  
  #get placebo results
  
  placebo.boot.cbw <- readRDS(paste0("results/placebo-boot-cbw-",o,cf,".rds"))[[1]] # Placebo t0 = 4
  
  testhat <- placebo.boot.cbw$t0 # test stat on placebo data
  test <- placebo.boot.cbw$t # boot statistics
  
  treat.status <- matrix(rownames(testhat), nrow=nrow(testhat), ncol=1)
  treat.status[rownames(testhat) %in% outcomes.cbw.placebo$eastern] <- "eastern"
  treat.status[rownames(testhat) %in% outcomes.cbw.placebo$swiss] <- "swiss"
  treat.status[rownames(testhat) %in% outcomes.cbw.placebo$control] <- "control"
  treat.status <- matrix(treat.status, dimnames=list(NULL, "status"))
  
  testhat.mean <-  aggregate(testhat, list(treat.status), mean)[-1]
  colnames(testhat.mean) <- colnames(outcomes.cbw.placebo$mask)
  
  testhat.eastern <-  rowMeans(testhat.mean[2,][,1:(which(colnames(outcomes.cbw.placebo$mask)=="20052"))]) # tau bar
  testhat.swiss <-  rowMeans(testhat.mean[3,][,1:(which(colnames(outcomes.cbw.placebo$mask)=="20054"))])
  
  test.eastern  <-apply(test, 1, function(x){
    x <- matrix(x, nrow=dim(testhat)[1], ncol=dim(testhat)[2], byrow=FALSE,
                dimnames=dimnames(testhat))
    test.mean <-  aggregate(x, list(treat.status), mean)[-1]
    colnames(test.mean) <- colnames(outcomes.cbw.placebo$mask)
    
    test.eastern <-  rowMeans(test.mean[2,][,1:(which(colnames(outcomes.cbw.placebo$mask)=="20052"))]) # tau bar
    
    return(test.eastern)
  })
  
  test.swiss  <-apply(test, 1, function(x){
    x <- matrix(x, nrow=dim(testhat)[1], ncol=dim(testhat)[2], byrow=FALSE,
                dimnames=dimnames(testhat))
    test.mean <-  aggregate(x, list(treat.status), mean)[-1]
    colnames(test.mean) <- colnames(outcomes.cbw.placebo$mask)
    
    test.swiss <-  rowMeans(test.mean[3,][,1:(which(colnames(outcomes.cbw.placebo$mask)=="20054"))])
    
    return(test.eastern)
  })
  
  p.val.eastern <- 1- ((1/length( test.eastern ) * sum( test.eastern  < testhat.eastern)))
  
  p.val.swiss <- 1- ((1/length( test.swiss ) * sum( test.swiss  < testhat.swiss)))
  
  return(list("eastern"=p.val.eastern,"swiss"=p.val.swiss))
})

names(outcome.pvals) <- outcome.vars
    
p.values <- data.frame("pvals"= c(unlist(lapply(outcome.pvals, '[[', 1)),unlist(lapply(outcome.pvals, '[[', 2))),
                       "Cluster"=c(rep("Eastern", length(outcome.vars)), rep("Swiss", length(outcome.vars))),
                       "Outcome"=rep(outcome.vars, 2))

#p.values.m <-melt(p.values,id.vars = c("Cluster","Outcome"))

#Plot
mc.placebo.plot <- ggplot(p.values, aes(x=pvals, y=Outcome)) + 
  geom_point(stat='identity', aes(col=factor(Cluster),shape=factor(Cluster)), size=4, alpha=0.9)  +
  labs(x="Placebo test p-values",
       y=" ") + 
  geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
  scale_x_continuous(breaks=c(0,0.05,0.1,0.25,0.5,0.75,1), 
                     labels=c("0","0.05","0.10","0.25","0.50","0.75","1")) +
  scale_y_discrete(labels=outcomes.labels, limits = levels(p.values$outcomes))+
  coord_flip() +
  scale_shape_manual(name="Cluster", values = c("Eastern" = 2,
                                                "Swiss" = 4),
                     labels=c("Eastern", "Swiss")) +
  scale_colour_manual(name="Cluster", values = c(  "Eastern" = wes_palette("Darjeeling1")[5],
                                                   "Swiss" = wes_palette("Darjeeling1")[4]),
                      labels=c("Eastern", "Swiss")) + theme(legend.key=element_blank(), legend.title=element_text(size=12)) + theme(plot.title = element_text(hjust = 0.5, size=12), plot.subtitle = element_text(hjust = 0.5, size=12)) + 
  theme(axis.text.y = element_text(size=12,margin = margin(t = 0, r = 0, b = 0, l = 20))) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8))

ggsave(filename = paste0("plots/placebo-pvals-cbw",cf,".png"),plot = mc.placebo.plot, scale=1.25)