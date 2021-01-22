getTrends <- function(L, plot=FALSE){
  svdL <- svd(L)
  svdd <- diag(svdL$d)
  svdd05 <- svdd^0.5
  V <- svdd05 %*% t(svdL$v)
  if(plot){
    plot(V[1,], type='l') #first latent trend
    plot(V[2,], type='l') #second latent trend
  }
  return(list("first"=V[1,],"second"=V[2,]))
}