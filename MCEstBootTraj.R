MCEstBootTraj <- function(impact,indices, eastern=NULL, swiss=NULL) {
  
  att <- impact[indices]

  if(!is.null(eastern)){
    return(mean(att[names(att) %in% eastern])) # avg over cluster
  }
  if(!is.null(swiss)){
    return(mean(att[names(att) %in% swiss])) # avg over cluster
  }
}