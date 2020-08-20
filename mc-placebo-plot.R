# Plot CIs from placebo tests
library(ggplot2)
library(wesanderson)
library(reshape2)

outcome.vars <- c("CBWbord","CBWbordEMPL","Thwusual")
outcomes.labels <- c("% working in border region",
                     "% working in border region,\n conditional on employment",
                     "Average total working hours")

T0 <- c("0.96875","0.8125","0.65625","0.5")

outcome.pvals <- lapply(T0, function(x){lapply(outcome.vars, function(o){lapply(c(0,1), function(i){

  #Analysis 1: ST vs AT (retrospective, X=CBW) 
  
  print(paste0("Estimates for outcome: ",o," Simultaneous: ",i))
  
  # placebo data
  outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  # Use post-treatment (all zeros)
  outcomes.cbw.placebo <- outcomes.cbw
  outcomes.cbw.placebo$mask <- outcomes.cbw$mask[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  outcomes.cbw.placebo$M <- outcomes.cbw$M[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  outcomes.cbw.placebo$W <- outcomes.cbw$W[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  
  #get placebo results
  
  placebo.boot.cbw <- readRDS(paste0("results/placebo-boot-cbw-",o,i,".rds"))
  
  testhat.eastern <- placebo.boot.cbw[[x]]$eastern$t0 # test stat on placebo data
  test.eastern <- placebo.boot.cbw[[x]]$eastern$t # boot statistics
  
  testhat.swiss <- placebo.boot.cbw[[x]]$swiss$t0 # test stat on placebo data
  test.swiss <- placebo.boot.cbw[[x]]$swiss$t # boot statistics
  
  p.val.eastern <- 1- ((1/length( test.eastern ) * sum( test.eastern  < testhat.eastern, na.rm = TRUE)))
  
  p.val.swiss <- 1- ((1/length( test.swiss ) * sum( test.swiss  < testhat.swiss, na.rm = TRUE)))
  
  print(paste0("Placebo T0: ", x))
  
  print(paste0("Pval Eastern: ", p.val.eastern))
  print(paste0("Pval Swiss: ", p.val.swiss))
})
})
})
