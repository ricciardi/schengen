packages <- c("dplyr", "readstata13", "glmnet", "caret", "ggplot2", "wesanderson", "reshape2", "zoo", "matrixStats", "tseries", "devtools", "Matrix", "tictoc", "MASS", "data.table", "reshape", "reshape2","scales")

install.packages(packages, repos = "http://cran.us.r-project.org")

install.packages("devtools")
library(devtools) 

install_github("jvpoulos/MCPanel")
install_github("Duane321/emfactor")

# doMPI
doMPI <- FALSE
if(doMPI){
  install.packages("Rmpi")
  install.packages("doMPI", dependencies=TRUE, repos = "http://cran.us.r-project.org")
}