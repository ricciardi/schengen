MCEstBootTraj <- function(impact,indices, t0.eastern=NULL, t0.swiss=NULL,eastern=NULL,swiss=NULL,start=1) {
  
  if(!is.null(t0.eastern)){
    att.eastern <- rowMeans(impact[indices,][,start:(t0.eastern-1)])# Schengen + FoM
    return(mean(att.eastern[names(att.eastern) %in% eastern])) # avg over cluster
  }
  if(!is.null(t0.swiss)){
    att.swiss <- rowMeans(impact[indices,][,start:(t0.swiss-1)]) # Schengen + FoM
    return(mean(att.swiss[names(att.swiss) %in% swiss])) # avg over cluster
  }
}