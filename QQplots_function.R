qqplots <- function(vector, radar_vector){
 
  # calculate the quantiles
  RD_q <- quantile(radar_vector, probs = seq(0.01, 0.99, 0.01), na.rm = T)
  k_q <- quantile(vector, probs = seq(0.01, 0.99, 0.01), na.rm = T)
  
  # calculate regression coefficient
  df <- data.frame(RD=as.vector(RD_q), k=as.vector(k_q))
  regr <- lm(k~RD, df)
  a <- coefficients(regr)[1]
  b <- coefficients(regr)[2]
  R2 <- 1-(sum((df$k - (a+b*df$RD))^2, na.rm=T)/sum((df$k - mean(df$k, na.rm=T))^2, na.rm=T))
  
  result <- data.frame(R2=R2, intercept = a, slope = b)

  return(result)
}
