MCEstBootTraj <- function(impact,indices, t0.eastern=NULL, t0.swiss=NULL) {
  
  if(!is.null(t0.eastern)){
    trajectory.eastern <- rowMeans(impact[,1:(t0.eastern-1)])[indices] # Schengen + FoM
    return(mean(att[names(att) %in% eastern])) # avg over cluster
  }
  if(!is.null(t0.swiss)){
    trajectory.swiss <- rowMeans(impact[,1:(t0.swiss-1)])[indices] # Schengen + FoM
    return(mean(att[names(att) %in% swiss])) # avg over cluster
  }
}