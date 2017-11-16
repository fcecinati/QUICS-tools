# Kriging with External Drift for Uncertain Data (KEDUD)
# kriging with external drift is used to merge radar and rain gauge data, conisdering spatially variable rain gauge errors:
# 
# Written by Francesca Cecinati
# University of Bristol
# Reference: 
#
# F. Cecinati, A. M. Moreno Ródenas, and M. A. Rico-Ramirez (2017), 
# "Integration of rain gauge errors in radar-rain gauge merging techniques," 
# in 10th World Congress on Water Resources and Environment, Athens, pp. 279 - 285.
#
# F. Cecinati, A. M. Moreno-Ródenas, M. A. Rico-Ramirez, M. ten Veldhuis, J. Langeveld (TBD), 
# "Considering rain gauge measurement uncertainty using Kriging for Uncertain Data", 
# Water Resources Management (in preparation).


KEDUD <- function(data, radar, newdata, variogr, model, errors){
  # data = SpatialPointsDataFrame with the RG measurements
  # radar = SpatialGridDataFrame with the radar measurements
  # newdata = SpatialGridDataFrame with the radar in the prediction points
  # variogr = variogram object
  # model = variogram model between 'Gau','Exp','Sph'
  # errors = vector of errors (same length of rain gauge data frame)
  
  # extract values from the variogram
  nugget <- variogr$psill[1]
  sill <- variogr$psill[2]
  range <- variogr$range[2]
  
  # rename
  names(data) <- c("rain")
  
  
  #covariance function
  if (model=='Exp'){
    cov_model <- function(d, n, s, r) {
      y <- d
      for (i in 1:length(d)){
        if (d[i]>0) {
          y[i] <- s -s*(1-exp(-3*d[i]/r))
        } else if (d[i]==0){
          y[i] = s+n
        }
      }
      return(y)
    }
  }
  if (model=='Gau'){
    cov_model <- function(d, n, s, r) {
      y <- d
      for (i in 1:length(d)){
        if (d[i]>0) {
          y[i] <- s -s*(1-exp(-3*(d[i]^2)/(r^2)))
        } else if (d[i]==0){
          y[i] = s+n
        }
      }
      return(y)
    }
  }
  if (model=='Sph'){
    cov_model <- function(d, n, s, r) {
      y <- d
      for (i in 1:length(d)){
        if (d[i]>0 & d[i]<=r) {
          y[i] <- s - s*(1.5*(d[i]/r)-0.5*((d[i]/r)^3))
        } else if (d[i]>r) {
          y[i] <- 0
        }else if (d[i]==0){
          y[i] = s+n
        }
      }
      return(y)
    }
  }

    
  # Distances of measurement points
  coord <- data@coords
  n <- dim(coord)[1]
  dis <- as.matrix(dist(coord,diag=FALSE, upper=FALSE))
  
  # extract radar data on RG locations
  extr <- over(data, radar)
  
  # Calculate covariance matrix
  C <- matrix(0,n,n)
  C <- apply(dis,MARGIN=c(1,2),FUN=cov_model, n=nugget, s=sill, r=range)
  C <- cbind(C,rep(1,n))
  C <- rbind(C,c(rep(1,n),0))
  C <- cbind(C,c(as.numeric(t(extr)), 0))
  C <- rbind(C,c(as.numeric(t(extr)),0,0))
  
  for (i in 1:n){
    C[i,i] <- C[i,i]+errors[i]^2
  }
  
  #Invert C
  Cinv <- solve(C)
  
  # Prepare matrices
  pred = matrix(0,newdata@grid@cells.dim[1],newdata@grid@cells.dim[2])
  var = matrix(0,newdata@grid@cells.dim[1],newdata@grid@cells.dim[2])
  predpoint = c(0,0)
  
  for (x in 1:newdata@grid@cells.dim[1]){
    for (y in 1:newdata@grid@cells.dim[2]) {
      
      predpoint <- coordinates(newdata[y,x])
      
      # Calculate distance at prediction point
      D<- matrix(1, 1, n+2)
      for (i in 1:n){
        D[i] <- ((coord[i,1]-predpoint[1])^2+(coord[i,2]-predpoint[2])^2)^(0.5)
      }
      D<-t(D)
      D<-apply(D,MARGIN=c(1,2),FUN=cov_model,n=nugget, s=sill, r=range)
      D[length(D)-1]<-1
      D[length(D)]<-newdata[y,x]@data[[1]]
      
      #Multiply inverse C and D to get the weights
      W <- Cinv%*%D
      
      # multiply the weights to the values
      w <- t(W[1:(length(W)-2)])
      pred[x,y] <- w%*%as.matrix(as.numeric(t(data@data)))
      
      # Variance calculation
      var[x,y] <- variogr[2]-t(W)%*%D
      #print(c(x,y))
    }
  }
  
  p <- raster(t(pred), xmn=newdata@bbox[1,1], ymn = newdata@bbox[2,1], xmx = newdata@bbox[1,2], ymx = newdata@bbox[2,2], crs = newdata@proj4string)
  v <- raster(t(var), xmn=newdata@bbox[1,1], ymn = newdata@bbox[2,1], xmx = newdata@bbox[1,2], ymx = newdata@bbox[2,2], crs = newdata@proj4string)
  p[is.na(newdata@data)] <- NA
  v[is.na(newdata@data)] <- NA
  p[p<0] <- 0
  v[v<0] <- 0
  res <- stack(p,v)
  results <- as(res, "SpatialGridDataFrame")
  
  return(results)
}