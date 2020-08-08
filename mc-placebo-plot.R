# Plot CIs from placebo tests
library(ggplot2)
library(wesanderson)
library(reshape2)

outcomes <- c("CBWbord","CBWbordEMPL","Thwusual")
outcomes.labels <- c("% working in border region",
                     "% working in border region,\n conditional on employment",
                     "Average total working hours")
cf<- ("")
T0 <- c("0.96875","0.8125","0.65625","0.5")

outcome.pvals <- lapply(T0, function(i){lapply(outcome.vars, function(o){
  print(i)
  print(o)
  
  #Analysis 1: ST vs AT (retrospective, X=CBW) 
  
  print(paste0("Estimates for Analysis 1, outcome:",o,cf))
  
  # placebo data
  outcomes.cbw <- readRDS(paste0("data/outcomes-cbw-",o,".rds"))
  
  # Use post-treatment (all zeros)
  outcomes.cbw.placebo <- outcomes.cbw
  outcomes.cbw.placebo$mask <- outcomes.cbw$mask[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  outcomes.cbw.placebo$M <- outcomes.cbw$M[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  outcomes.cbw.placebo$W <- outcomes.cbw$W[,which(colnames(outcomes.cbw$mask)=="20111"):ncol(outcomes.cbw$mask)]
  
  #get placebo results
  
  placebo.boot.cbw <- readRDS(paste0("results/placebo-boot-cbw-",o,"1.rds"))
  
  testhat.eastern <- placebo.boot.cbw[[i]]$eastern$t0 # test stat on placebo data
  test.eastern <- placebo.boot.cbw[[i]]$eastern$t # boot statistics
  
  testhat.swiss <- placebo.boot.cbw[[i]]$swiss$t0 # test stat on placebo data
  test.swiss <- placebo.boot.cbw[[i]]$swiss$t # boot statistics
  
  p.val.eastern <- 1- ((1/length( test.eastern ) * sum( test.eastern  < testhat.eastern, na.rm = TRUE)))
  
  p.val.swiss <- 1- ((1/length( test.swiss ) * sum( test.swiss  < testhat.swiss, na.rm = TRUE)))
  
  return(list("eastern"=p.val.eastern,"swiss"=p.val.swiss))
})
})

outcome.pvals