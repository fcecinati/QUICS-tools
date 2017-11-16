gaussianity <- function(vector){

  ###################################### load libraries  ####################################
  library(e1071) 
  source("Negentropy.R")
  ###########################################################################################
  
  kur <- kurtosis(as.vector(t(vector)), na.rm=T, type=2)
  skw <- skewness(as.vector(t(vector)), na.rm=T, type=2)
  neg <- negentropy(vector)
  
  
  result <- c(kurtosis=kur, skweness=skw, negentropy=neg)
  return(result)
}
