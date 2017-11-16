RG_validation <- function(ob, kd){
  
  df <- data.frame(ob=ob, kd=kd)
  df <- df[complete.cases(df),]
  ob <- as.numeric(df$ob)
  kd <- as.numeric(df$kd)
  
  MRTE <- mean((sqrt(kd)-sqrt(ob))^2)
  BIAS <- mean(abs(kd-ob))
  
  ob <- round(ob*5)/5
  kd <- round(kd*5)/5
  
  RG_bin <- ob; RG_bin[RG_bin>0]<-1
  ked_bin <- kd; ked_bin[ked_bin>0] <- 1
  A <- sum(RG_bin==1 & ked_bin==1)
  B <- sum(RG_bin==0 & ked_bin==1)
  C <- sum(RG_bin==1 & ked_bin==0)
  D <- sum(RG_bin==0 & ked_bin==0)
  
  HK <- ((A*D)-(B*C))/((A+C)*(B+D))
  
  result <- c(MRTE, BIAS, HK)
  return(result)
  
}