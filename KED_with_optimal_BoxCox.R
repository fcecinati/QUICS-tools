# Kriging with External Drift with Optimised time-variant Box-Cox Transformation 
# applied to Radar and Rain Gauge Merging
# 
# Written by Francesca Cecinati
# University of Bristol
# Reference: 
#
#	F. Cecinati, O. Wani, and M. A. Rico-Ramirez (2017), 
# "Comparing approaches to deal with non-Gaussianity of rainfall data in Kriging-based radar-gauge rainfall merging" 
# Water Resour. Res., (Accepted), DOI: 10.1002/2016WR020330.


#   !  !  !         REQUIRES:           !  !  !         
# Residual_gaussianity_function.R, QQplots_function.R, RG_validation_function.R, Negentropy.R


#########################################################################################################################
rm(list = ls())  # clean memory
gc()
#######################################           Input              ###################################################

# Dataset files:
RGfile <- "X:/UK/P3/Dataset/AfterReview/EA"                       # Table T x N where each line is a time step and each column a location
RGfilecoord <- "X:/UK/P3/Dataset/AfterReview/EA_coord"            # Table N x 2 with columns corresponding to rain gauge coordinates X and Y
RDfile <- "X:/UK/P3/Dataset/AfterReview/RD"                       # Table T x M where each line is a time step and each column a grid point
RDfilecoord <- "X:/UK/P3/Dataset/AfterReview/RD_coord"            # Table M x 2 with columns corresponding to radar grid coordinates X and Y
RDfiledates <- "X:/UK/P3/Dataset/AfterReview/RD_dates"            # Table T x 6 where each line is a time step and columns are Y m D H M S

# Variogram destination folder
variogramfolder <- "X:/UK/P3/Variogram/AfterReview/"

# Projection parameters in proj4 format
myproj <- "+proj=stere +lat_0=90 +lon_0=0.0 +lat_ts=60.0 +a=6378.137 +b=6356.752 +x_0=0 +y_0=0"

# Result folder
Results_ked <- paste0("X:/UK/P3/Results/AfterReview/BC_opt/")

#######################################            Functions        #####################################################

# Semivariogram model
v_model <- function(d,n,s,r){
  y <- d
  for (i in 1:length(d)){
    if (d[i]>0) {
      y[i] <- n + (s)*(1-exp(-3*d[i]/r))
    } else if (d[i]==0){
      y[i] = 0
    }
  }
  return(y)
}

# V_model least squares:
v_ls <- function(param,d,variog, weights){
  n <- param[1]
  s <- param[2]
  r <- param[3]
  mod <- v_model(d, n, s, r)
  obs <- variog
  ls <- weighted.mean((obs-mod)^2, weights, na.rm=T)
  return(ls)
}

# Box-Cox function
bc <- function(y,l){
  if (l==0){
    bc <- log(y)
  } else {
    bc <- (y^l - 1)/l
  }
  return(bc)
}

# Box-Cox back function
bcb <- function(y,l){
  if (l==0){
    bcb <- exp(y)
  } else {
    bcb <- (l*y+1)^(1/l)
  }
  return(bcb)
}

###################################### load libraries  ####################################
library(sp)
library(maptools)
library(gstat)
library(spacetime)
library(raster)
library(minpack.lm)
library(lubridate)
library(plyr)
library(gdata)
library(ggplot2)
library(plotKML)
library(gridExtra)
library(RSAGA)
library(abind)
library(hydromad)
source("Residual_gaussianity_function.R")
source("QQplots_function.R")
source("RG_validation_function.R")
source("Negentropy.R")
###########################################################################################
# Read files

# Rain gauge files
load(RGfile)
load(RGfilecoord)

# read radar
load(RDfile)
load(RDfilecoord)
load(RDfiledates)

# Dimensions
t_steps <- dim(RGtable)[1]
n_RG <- dim(RGtable)[2]
n_RD <- dim(RDtable)[2]

# Prepare output
if (file.exists(paste0(Results_ked, "KED.txt"))){
  file.remove(paste0(Results_ked, "KED.txt"))
}

# Preapre the data
var_KED_BC <- data.frame(matrix(NA, dim(RGtable)[1], 3)); names(var_KED_BC) <- c("nugget", "psill", "range")

test1 <- data.frame(matrix(NA, dim(RGtable)[1], 4)); names(test1) <- c("kurtosis", "skewness", "negentropy", "kolmogorov-smirnov")
test2 <- data.frame(matrix(NA, dim(RGtable)[1], 3)); names(test2) <- c("R2", "intercept", "slope")
test3 <- data.frame(matrix(NA, dim(RGtable)[1], 3)); names(test3) <- c("MRTE", "bias", "HK")

RG_valid <- matrix(NA, dim(RGtable)[1], dim(RGtable)[2])
residuals <- matrix(NA, dim(RGtable)[1], dim(RGtable)[2])

lambda=rep(NA, t_steps)

# Loop on the available time steps
for (t in 1:t_steps){

  # prepare one rain gauge tstep
  RG <- RGtable[t,]
  date <- row.names(RG); date <- strptime(date, format="X%Y.%m.%d.%H.%M.%S", tz="UTC")
  date_index <- which((RDdates[,1]==year(date))&(RDdates[,2]==month(date))&(RDdates[,3]==day(date))&(RDdates[,4]==(hour(date)-1)))
  RGdates <- c(year(date), month(date), day(date), hour(date), 0, 0)


  if ((mean(unlist(RGtable[t,]), na.rm=T)>0)&(length(date_index)>0)){

    # RG as data frame
    RG <- data.frame(t(RG), RGcoord)
    names(RG) <- c("rain", "x", "y")
    coordinates(RG) <- ~x+y
    proj4string(RG) <- myproj

    # radar as data frame
    RD <- RDtable[date_index,]
    RD <- data.frame(t(RD), RDcoord)
    names(RD) <- c("radar", "x", "y")
    coordinates(RD) <- ~x+y
    proj4string(RD) <- myproj
    RD <- as(RD, "SpatialPixelsDataFrame")
    
    # Extract radar at rain gauge locations
    RG$radar <- over(RG, RD)$radar
    
    # Perform merging only when sufficient rain is recorded (adjust parameters)
    if ((mean(unlist(RG$rain), na.rm=T)>0) & (mean(unlist(RG$radar), na.rm=T)>0.001) & sum(RD@data$radar>0, na.rm=T)>400){

      lambdas <- seq(from=0.20, to=1.50, by=0.01 )
      ng <- rep(NA, length(lambdas))
      
      # identify a time window where at least 500 points are available
      nobs <- 0
      window <- 0
      while (nobs<500) {
        window <- window + 1
        nobs <- sum(is.na(RGtable[max((t-window),1):min((t+window), t_steps),])==F)
      }
      RG_sub <- RGtable[max((t-window),1):min((t+window), t_steps),]
      RG_vector <- NA
      RD_vector <- NA
      
      # collact all teh available data in the time window
      for (ii in 1:dim(RG_sub)[1]){
        
        # preapre rain gauges
        RGs <- data.frame(t(RG_sub[ii,]), RGcoord)
        names(RGs) <- c("rain", "x", "y")
        id_RG <- complete.cases(RGs)
        RGs <- RGs[complete.cases(RGs),]
        coordinates(RGs) <- ~x+y
        proj4string(RGs) <- myproj
        
        # identify corresponding radar
        dt <- row.names(RG_sub[ii,]); dt <- strptime(dt, format="X%Y.%m.%d.%H.%M.%S", tz="UTC")
        dtind <- which((RDdates[,1]==year(dt))&(RDdates[,2]==month(dt))&(RDdates[,3]==day(dt))&(RDdates[,4]==(hour(dt)-1)))
        
        # if radar is available for the time step
        if (length(dtind)>0){
          # radar as data frame
          RDs <- RDtable[dtind,]
          RDs <- data.frame(t(RDs), RDcoord)
          names(RDs) <- c("radar", "x", "y")
          coordinates(RDs) <- ~x+y
          proj4string(RDs) <- myproj
          RDs <- as(RDs, "SpatialPixelsDataFrame")
          RGs$radar <- over(RGs, RDs)$radar
          
          # Line up all the available data over the different time steps
          RG_vector <- c(RG_vector, RGs$rain)
          RD_vector <- c(RD_vector, RGs$radar)
        }
      }

      RG_vector <- RG_vector[is.na(RG_vector)==F]
      RD_vector <- RD_vector[is.na(RD_vector)==F]

      # Optimization of lambda
      for (i in 1:length(lambdas)){

        # calculation of residuals with transformation
        rsd <- bc(RG_vector, lambdas[i]) - bc(RD_vector, lambdas[i])
        ng[i] <- negentropy(rsd)
      }
      # optional plot
      # plot(lambdas, ng)

      # select the lambda that minimises the negentropy
      lambda[t]=lambdas[ng==min(ng)]

      # Return to numeric vectors
      RG <- RGtable[t,]
      RD <- RDtable[date_index,]

      # Apply Box-Cox transformation
      l=lambda[t]
      if (is.na(l)==T){
        l=0.5
      }
      RG <- bc(RG, l)
      RD <- bc(RD, l)

      #############################################################################
      ############                Calculate variogram               ###############
      #############################################################################

      # data frame
      RDv <- data.frame(t(RD), RDcoord)
      RDv[RDv==bc(0,l)] <- NA
      RDv <- RDv[complete.cases(RDv),]
      names(RDv) <- c("rain", "x", "y")

      # Use a subset of the radar > 400 and < 1000 wet points to calculate the variogram
      if (dim(RDv)[1]>1000){
        n_RD_sub <- 1000
      } else {
        n_RD_sub <- dim(RDv)[1]
      }

      #use only a subset of the radar
      sub_index <- round(runif(n_RD_sub, min=1, max=dim(RDv)[1]))
      RD_sub <- RDv[sub_index,]

      # Calculate distances
      dis <- matrix(0, n_RD_sub, n_RD_sub)
      disx <- dist(RD_sub$x)
      disy <- dist(RD_sub$y)
      dis2 <- sqrt(disx^2 + disy^2)
      dis[lower.tri(dis)] <- dis2

      # Empirical rainfall Variogram
      maxd <- ceiling(max(dis))
      d <- seq(from=1000, to=maxd, by=1000)
      variog <- d
      num <- d
      c=0
      for (i in d){
        c=c+1
        indices <- which((dis>i-1000)&(dis<=i), arr.ind=TRUE)

        if (dim(indices)[1]>0){
          mat <- (RD_sub[indices[,1], "rain"]-RD_sub[indices[,2], "rain"])^2
          variog[c] <- mean(mat, na.rm=TRUE)/2
          num[c] <- sum(is.finite(mat))
        } else {
          variog[c] <- NA
          num[c] <- NA
        }
      }
      d <- d[1:floor(length(d)/2)]
      variog <- variog[1:floor(length(d))]
      num <- num[1:floor(length(d))]
      weights <- num/(d^2)  #weight for the fitting

      # in case of outliers
      if (mean(variog, na.rm=T)>1000){
        nugget=0
        range=40000
        sill=1
      } else {
        # Finding the optimum variogram parameters with least square minimizaton with SCEoptim (shuffle complex evolution)
        st_guess <- c(n=0.1, s=mean(variog, na.rm=T), r=maxd/8)
        lbound <- c(n=0, s=0.000001, r=0)
        ubound <- c(n=1, s=200, r=500000)
        ctr <- list(reltol=10E-7, maxit = 20, trace = F, ncomplex=6)
        optPars <- SCEoptim(
          FUN=v_ls,
          par=st_guess, lower=lbound, upper=ubound,
          control=ctr,
          d=d, variog=variog, weights=weights
        )

        nugget <- optPars$par[1]
        sill <- optPars$par[2]
        range <- optPars$par[3]

      }

      #Save the OK variogram
      vgm_ok <- vgm(psill=sill, nugget=nugget, range=range, model="Exp")

      ##### Calculate residuals

      # radar as data frame
      RD <- data.frame(t(RD), RDcoord)
      names(RD) <- c("radar", "x", "y")
      coordinates(RD) <- ~x+y
      proj4string(RD) <- myproj
      RD <- as(RD, "SpatialPixelsDataFrame")
      RD <- as(raster(RD), "SpatialGridDataFrame")

      # RG as data frame
      RG <- data.frame(t(RG), RGcoord)
      names(RG) <- c("rain", "x", "y")
      id_RG <- complete.cases(RG)
      RG <- RG[complete.cases(RG),]
      coordinates(RG) <- ~x+y
      proj4string(RG) <- myproj

      # Ordinary Kriging
      RG_ok <- krige(rain~1, RG, newdata=RD, model=vgm_ok, na.action = na.pass)
      names(RG_ok) <- c("pred", "var")
      RG_ok$pred[RG_ok$pred<bc(0,l)] <- bc(0,l)

      # residuals
      df <- data.frame(RD=RD@data$radar, OK=RG_ok$pred)
      resid_lm <- lm(OK~RD, df)
      resid <- RD
      resid@data$rain <- RG_ok@data$pred - coef(resid_lm)[2]*RD@data$radar+coef(resid_lm)[1]

      #use only a subset of the residuals
      res_sub <- data.frame(resid)
      res_sub <- res_sub[sub_index,]

      # Empirical residual
      c=0
      for (i in d){
        c=c+1
        indices <- which((dis>i-1000)&(dis<=i), arr.ind=TRUE)
        if (dim(indices)[1]>0){
          mat <- (res_sub[indices[,1], "rain"]-res_sub[indices[,2], "rain"])^2
          variog[c] <- mean(mat, na.rm=TRUE)/2
          num[c] <- sum(is.finite(mat))
        } else {
          variog[c] <- NA
          num[c] <- NA
        }
      }
      d <- d[1:floor(length(d)/2)]
      variog <- variog[1:floor(length(d))]
      num <- num[1:floor(length(d))]

      # Rainfall variogram fitting
      weights <- num/(d^2)  #weight for the fitting

      if (mean(variog, na.rm=T)>1000){
        nugget=0
        range=40000
        sill=1
      } else {
        # Finding the optimum variogram parameters with least square minimizaton with SCEoptim (shuffle complex evolution)
        st_guess <- c(n=0.1, s=mean(variog, na.rm=T), r=maxd/8)
        lbound <- c(n=0, s=0.000001, r=0)
        ubound <- c(n=1, s=200, r=500000)
        ctr <- list(reltol=10E-7, maxit = 20, trace = F, ncomplex=6)
        optPars <- SCEoptim(
          FUN=v_ls,
          par=st_guess, lower=lbound, upper=ubound,
          control=ctr,
          d=d, variog=variog, weights=weights
        )

        nugget <- optPars$par[1]
        sill <- optPars$par[2]
        range <- optPars$par[3]

      }

      #Save the residuals variogram
      var_KED_BC[t,] <- c(nugget,sill,range)
      vgm_ked <- vgm(psill=sill, nugget=nugget, range=range, model="Exp")


      #############################################################################
      ############      Test 1 : gaussianity of the residuals       ###############
      #############################################################################
      
      test1[t,] <- as.numeric(gaussianity(resid@data$rain))

      #############################################################################
      ############                Perform merging                   ###############
      #############################################################################

      # attach the radar to the rain gauge data frame
      RG$radar <- over(RG, RD)$radar

      # Calculate residuals
      residuals[t,is.na(RGtable[t,])==F] <- RG$rain - RG$radar

      # KED
      ked <- krige(rain~radar, RG, newdata=RD, vgm_ked, na.action = na.pass)
      names(ked) <- c("pred", "var")

      k <- matrix(ked$pred, 201, 201, byrow=T)
      k <- k[dim(k)[1]:1,]
      k <- as.numeric(k)

      va <- matrix(ked$var, 201, 201, byrow=T)
      va <- va[dim(va)[1]:1,]
      va <- as.numeric(va)

      # Back-transform
      quant <- matrix(NA, length(k), 99)
      bquant <- matrix(NA, length(k), 99)
      for (q in 1:99){
        quant[,q] <- qnorm(q/100, mean=k, sd=sqrt(va))
        bquant[,q] <- bcb(quant[,q], l)
        bquant[is.na(bquant)==T] <- 0 # when the value falls below a certain threshold, there is a negative number in the back transform that is elevated to a certain power. corresponds to 0
      }
      k_bcb <- apply(bquant, 1, mean)
      k_bcb[k_bcb<0] <- 0
      k_bcb <- round(t(as.numeric(c(RGdates, k_bcb))), digits = 4)

      write.table(k_bcb, file=paste0(Results_ked, "KED.txt"), append=T, quote = F, sep="\t", row.names=F, col.names=F)

      ###############################################################################
      ############   Test 2 : reproduction of original distribution   ###############
      ###############################################################################

      vector <- k_bcb[7:length(k_bcb)]
      radar_vector <- RDtable[date_index,]
      test2[t,] <- as.numeric(qqplots(vector, radar_vector))

      ###############################################################################
      ############          Cross-validation for KED                  ###############
      ###############################################################################

      # loop over the available RG
      nRG <- dim(RG)[1]
      prediction <- rep(NA, nRG)
      for (j in 1:nRG){

        # divide the dataset
        RGdf <- data.frame(RG)
        RGval <- RGdf[j,]
        RGpred <- RGdf[-j,]

        coordinates(RGval) <- ~x+y
        proj4string(RGval) <- myproj
        coordinates(RGpred) <- ~x+y
        proj4string(RGpred) <- myproj

        # KED only on RGval
        ked_val <- krige(rain~radar, RGpred, newdata=RGval, vgm_ked, na.action = na.pass)
        names(ked_val) <- c("pred", "var")

        p <- ked_val$pred
        v <- ked_val$var

        # Back-transform
        quant <- rep(NA, 99)
        bquant <- rep(NA, 99)
        for (q in 1:99){
          quant[q] <- qnorm(q/100, mean=p, sd=sqrt(v))
          bquant[q] <- bcb(quant[q], l)
          bquant[is.na(bquant)==T] <- 0 # when the value falls below a certain threshold, there is a negative number in the back transform that is elevated to a certain power. corresponds to 0
        }
        p_bcb <- mean(bquant)

        prediction[j] <- p_bcb
      }

      prediction[prediction<0] <- 0
      RG_valid[t,id_RG] <- prediction


      ###############################################################################
      ############            Test 3 : validation for KED             ###############
      ###############################################################################
      #prediction <- round(prediction*5)/5
      ob <- as.numeric(RGtable[t,])
      kd <- as.numeric(RG_valid[t,])
      test3[t,] <- as.numeric(RG_validation(ob, kd))

    }else {
      # just write NAs
      write.table(t(as.numeric(c(RGdates, rep(NA, dim(RDtable)[2])))), file=paste0(Results_ked, "KED.txt"), append=T, quote = F, sep="\t", row.names=F, col.names=F)
    }
  } else {
    # just write NAs
    write.table(t(as.numeric(c(RGdates, rep(NA, dim(RDtable)[2])))), file=paste0(Results_ked, "KED.txt"), append=T, quote = F, sep="\t", row.names=F, col.names=F)
  }
}

# Plot lambda
lambda <- data.frame(lambda)
names(lambda) <- c("l")

tiff(filename=paste0(Results_ked, "Lambda_hist.tiff"), height=2000, width=2500, res=300)
p <- ggplot(lambda, aes(l))+
  geom_histogram(colour="darkred", fill="orangered", binwidth=0.02)+
  theme(axis.text=element_text(size=18, colour="azure4"),
        axis.title=element_text(size=18),
        panel.background=element_rect(fill="white"),
        panel.grid.major=element_line(colour="slategray2"),
        panel.grid.minor=element_line(colour="slategray2"),
        plot.title=element_text(size=24, face="bold", hjust=0.5))+
  ggtitle("Optimal l")+
  xlab("l value")+
  ylab("Frequency")
print(p)
dev.off()

save(lambda, file=paste0(Results_ked, "lambda"))
save(var_KED_BC, file=paste0(variogramfolder, paste0("KED_BCopt")))
save(test1, file=paste0(Results_ked, "test1"))
save(residuals, file=paste0(Results_ked, "residuals"))
save(test2, file=paste0(Results_ked, "test2"))
save(test3, file=paste0(Results_ked, "test3"))
save(RG_valid, file=paste0(Results_ked, "RG_valid"))
save(window_count, file=paste0(Results_ked, "window"))
