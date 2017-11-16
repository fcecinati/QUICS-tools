################################################################################################################
#      This function calculates an approximation of negentropy as in Hyvarinen, 1998, for a given vector      #
################################################################################################################
negentropy <- function(vector){

  # Function to use for approximated negentropy
  G_fun <- function(y){
    G <- -exp(-(y^2)/2)
    return(G)
  }

  
  # if it is not a vector, convert into a vector
  vector = as.vector(t(vector))
  
  # Calculate variable stats:
  mu <- mean(vector, na.rm=T)
  var <- var(vector, na.rm=T)
  std <- sqrt(var)
  v_length <- length(vector)
  
  # Standardise the variable
  vector_st <- (vector-mu)/std
  
  # Mean of the function G applied to a large normal sample
  mG_rnd <- -0.7071
  
  # Calculate approximation of negentropy
  G_VEC <- G_fun(vector_st)
  
  J <- (mean(G_VEC, na.rm=T)-mG_rnd)^2
  
  return(J)
}
