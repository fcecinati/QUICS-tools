# Geostatistical Radar Ensemble Generator (GREG)
# 
# Written by Francesca Cecinati
# University of Bristol
# Reference: DOI:10.1016/j.jhydrol.2017.02.053

#########################################################################################################################
rm(list = ls())  # clean memory
gc()

######################################        load libraries and functions          #####################################
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
library(RSAGA)
library(abind)
library(foreach)
library(doMC)
registerDoMC(4)  #change to your number of CPU cores  


# Variogram model (exponential)
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
########################################           Input              ###################################################
#Starting date:
starting = '2007/10/01 00:00'

# Results path
results_path = '/home/fc14509/Paper1/Ensembles/'

# Results path for validation-extracted ensembles
val_path <- '/home/fc14509/Paper1/Extracted ensembles2/Rain Gauge validation/'
hydro_path <- '/home/fc14509/Paper1/Extracted ensembles2/Catchment Ensembles/'
variog_path <- '/home/fc14509/Paper1/Extracted ensembles2/Variogram_parameters/'

# Path to the complete radar data in NIMROD format
radar_folder = '/home/fc14509/UK/Data/NIMROD_hourly/DAT/'  ##### NOTE: store the radar data in a folder inside this path named with the year e.g. "2007"

# Temporary folder to decompress the files
tmpfolder = '/home/fc14509/UK/Temp/'

# Path to the folder where the rain gauges, the radar data extracted at rain gauge location, and the coordinate .dat files are:
data_path = '/home/fc14509/data_north_england/'

# Name of the coordinate files of the observation points (two columns: x, Y):
coordinates_file = 'gaugecoordinates.dat'

# Name of the rain gauges .dat file (first three columns are (Y m D H M S), then each line are the data for all the grid points as per coordinate file:
rain_gauge_file = 'gauge_1h_2007-2011.dat'

# Name of the corresponding .dat radar file extracted at rain gauge location in the same format:
corresponding_radar_file = 'radar_1h_2007-2010.dat'

# Three Catchments shape files:
Rawtheyfile <- "/home/fc14509/PDM/Catchments/Rawthey/Rawthey.shp"
Ribblefile <- "/home/fc14509/PDM/Catchments/Ribble/Ribble.shp"
Lunefile <- "/home/fc14509/PDM/Catchments/Lune/Lune.shp"

# validation rain gauges:
validation = c(2,4,10,16,17,25,48,62,65,68,71,106,115,158,161,189)

# Projection system in proj4string format
myproj <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs"


# Parameters:

# How many ensembles do you want to generate
n_ensembles = 100

# How many hours is your simulation long?
t_steps = 8784 

# what are the easting and nothing in meters (top left corner)?
easting = 320000
northing = 550000

# What is the number of pixels for the square side?
side = 180000

#####################################                 Start script                #######################################
  
#Results folder
if (!dir.exists(results_path)){
  dir.create(results_path)
}

# Read rain gauge data
coord <- read.table(paste0(data_path, coordinates_file)); names(coord) <- c("X", "Y")
G <- read.table(paste0(data_path,rain_gauge_file)) 
R <- read.table(paste0(data_path,corresponding_radar_file))

# Separate the date part
G_start <- which(G[,1]==year(starting) & G[,2]==month(starting) & G[,3]==day(starting) & G[,4]==hour(starting))
G_end <- G_start + t_steps - 1
dates <- G[G_start:G_end, 1:6]
G <- G[G_start:G_end,7:dim(G)[2]]

# Separate the validation rain gauges
G_val <- G[,validation]; coord_val <- coord[validation,]
names(G_val) <- rownames(coord_val)

# Read radar data at rain gauge location
R_start <- which(R[,1]==year(starting) & R[,2]==month(starting) & R[,3]==day(starting) & R[,4]==hour(starting))
R_end <- R_start + t_steps - 1
R <- R[R_start:R_end,7:dim(R)[2]]
R_val <- R[,validation]; names(R_val) <- rownames(coord_val)

#Use of NA for negative values and to ignore days of no rain
G0 <- G
G[G<0.1] <- NA
R[R<0.1] <- NA
G0[G0<0] <- NA             # The version with zeroes need to be used in the variance correction
G_val[G_val<0] <- NA
R_val[R_val<0] <- NA

save(G_val, R_val, coord_val, file=paste0(val_path,'G_val'))
  
#Final size of the data matrix
x=dim(G)[2]

# Calculation of error statistics
er <- 10*(log10(R/G))
er <- matrix(unlist(er),t_steps,x)
mean_e <- mean(er, na.rm=T)
std_e <- sd(er,na.rm=T)
  
save(G, R, coord, er, file=paste0(results_path,'Data'))
save(dates, file=paste0(results_path, 'Dates'))
  
# Calculate distances
dist <- matrix(0, x, x)
for (i in 1:x){
  for(j in 1:x){
    dist[i,j] <- sqrt((coord[i,1]-coord[j,1])^2 + (coord[i,2]-coord[j,2])^2)
  }
}
dist <- dist*(lower.tri(dist, diag = FALSE))

# General Empirical Variogram on all available dates
maxd <- ceiling(max(dist))

d <- seq(from=1000, to=maxd, by=1000) 
variog <- d
num <- d
c=0
for (i in d){
  c=c+1
  indices <- which((dist>i-1000)&(dist<=i), arr.ind=TRUE)
  
  if (dim(indices)[1]>0){
    mat <- (er[,indices[,1]]-er[,indices[,2]])^2
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

# Variogram fitting
weights <- num/(d^2)  #weight for the fitting
myfit <- nlsLM(variog ~ v_model(d,n,s,r), data=data.frame(d,variog), weights=weights, na.action=na.omit, start=c(n=0,s=10, r=40000), lower = c(n=0, s=0, r=0), upper = c(n=200, s=200, r=10000000))
nugget <- coefficients(myfit)[1]
sill <- coefficients(myfit)[2]
range <- coefficients(myfit)[3]

variogram <- c(nugget,sill,range)
save(variogram, file=paste0(results_path, "General_Variogram"))

# Time variant part

# Variables pre-allocation
time <- 0
tr <- 0
class(tr)="try-error"
ens_val <- array(data=NA, dim=c(dim(G_val)[2],100))
ens_hydro <- array(data=NA, dim=c(3,100))

# Read shape files
Lune <- readShapePoly(Lunefile)
Rawthey <- readShapePoly(Rawtheyfile)
Ribble <- readShapePoly(Ribblefile)
Lune@proj4string@projargs <- myproj
Rawthey@proj4string@projargs <- myproj
Ribble@proj4string@projargs <- myproj

# Time loop in parallel
foreach (t=1:t_steps) %dopar% {    
  
  # variable reset
  n_hours=0
  sill_spec=0
  range_spec=0
  nugget_spec=0
  stdev=0
  mu=0
  availableobservations=0
  missing=0
  count_bignugget=0
  
  # use a try check for variogram fitting
  while (class(tr)=="try-error" & n_hours<12){
    
    # increase the number of considered time steps for the variogram calculation
    n_hours <- n_hours+3
    
    # check if we can subtract 3 hours, in order to distinguish the starting time steps
    if (t-n_hours<=0){
      t_start=1
    } else{
      t_start=t-n_hours+1
    }
    
    errors = er[t_start:t,]; errors <- matrix(errors, (t-t_start+1), dim(er)[2])
    maxd <- ceiling(max(dist))
    
    # Specific empirical Variograms
    d <- seq(from=1000, to=maxd, by=1000) 
    variog <- d
    num <- d
    c=0
    for (i in d){
      c=c+1
      indices <- which((dist>i-1000)&(dist<=i), arr.ind=TRUE)
      
      if (dim(indices)[1]>0){
        mat <- (errors[,indices[,1]]-errors[,indices[,2]])^2
        variog[c] <- mean(mat, na.rm=TRUE)/2
        num[c] <- sum(is.finite(mat))
      } else {
        variog[c] <- NA
        num[c] <- NA
      }
    }
    
    # reduce the considered variogram points to half of the maximum distance
    d <- d[1:floor(length(d)/2)]
    variog <- variog[1:floor(length(d))]
    num <- num[1:floor(length(d))]
    
    # Variogram fitting
    weights <- num/(d^2)  #weight for the fitting
    availableobservations = sum(num, na.rm=t)
    tr <- try(myfit <- nlsLM(variog ~ v_model(d,n,s,r), data=data.frame(d,variog), weights=weights, na.action=na.omit, start=c(n=0,s=mean(variog, na.rm=T), r=40000), lower = c(n=0, s=0, r=0)), silent=T)
    nugget_spec <- coefficients(myfit)[1]
    sill_spec <- coefficients(myfit)[2]
    range_spec <- coefficients(myfit)[3]
    stdev <- sd(errors, na.rm=T)
    mu <- mean(errors, na.rm=T)
    if (availableobservations<100 | nugget_spec>=sill_spec){   # even if it fits, if the minimum requirements are not satisfied it consider it failed
      class(tr)="try-error"
      if (nugget_spec>=sill_spec){
        count_bignugget<- 1
      }
    }
  }
  
  # either if it fits or jumps out after the 12 hour try
  
  if (class(tr)=="try-error"){ # case in which it fails to fit we use the general variogram
    sill_spec <- sill
    range_spec <- range
    nugget_spec <- nugget
    stdev <- std_e
    mu <- mean_e
    n_hours <- 15
  } else if (availableobservations<100 | nugget_spec>=sill_spec){ # case in which it does not satisfy the requirements we use the general variogram
    sill_spec <- sill
    range_spec <- range
    nugget_spec <- nugget
    stdev <- std_e
    mu <- mean_e
    n_hours <- 15
  } else {  #case in which the fitting is successful
    nugget_spec <- coefficients(myfit)[1]
    sill_spec <- coefficients(myfit)[2]
    range_spec <- coefficients(myfit)[3]
    stdev <- sd(errors, na.rm=T)
    mu <- mean(errors, na.rm=T)
  }

  class(tr)="try-error" # for the next time step, to enter in the while loop
  
  # open the NIMROD radar file (check the file naming)
  radarname <- paste0(toString(dates[t,1]), '/',toString(dates[t,1]),sprintf("%02d", dates[t,2]),sprintf("%02d", dates[t,3]),sprintf("%02d", dates[t,4]),sprintf("%02d", dates[t,5]))
  
  if (!file.exists(paste0(radar_folder, radarname))){ # if the radar file does not exist
    missing <- 1
  } else{
    variog <- vgm(psill=sill_spec, range=range_spec, nugget=nugget_spec, model='Exp')
    
    # Prepare radar data
    r <- read.table(paste0(radar_folder,radarname), sep=';'); r <- matrix(unlist(r), dim(r)[1], dim(r)[2])
    param <- read.table(paste0(radar_folder,'parameters'), sep=';'); param <- unlist(param)
    
    # invert the vertical dimension
    r <- r[1000:1,] 
    
    # prepare the spatial grid data frame
    r <- raster(r, xmn=(param[3]), xmx=(param[3] + param[7]*param[5]),
                ymn=(param[2] - param[4]*param[6]), ymx=(param[2]), 
                crs=myproj)
    e  <- extent(easting, easting+side, northing-side, northing) 
    r <- crop(r,e)
    r <- as(r, "SpatialGridDataFrame")
    
    # Error dataframe
    err <- errors
    err <- apply(err, 2, mean, na.rm=T)
    err <- data.frame(err,coord)
    err <- err[complete.cases(err),]

    
    # Check if there are observations to make a conditional simulation on
    if (dim(err)[1]==0){ # in case there are no observations, the simulation is unconditional
      
      #define prediction grid
      datagrid <- r
      datagrid$layer <- datagrid$layer*0
      datagrid@proj4string@projargs <- myproj
      
      # unconditional simulation
      g.obj <- gstat(formula=z~1, locations = ~x+y, dummy=T, beta=mu, model=variog, nmax=24)
      sim <- predict(g.obj, datagrid, nsim = 100)
      
      # Generate the ensembles
      log_ens <- sim
      ens_nc <- sim
      for (i in 1:100){
        log_ens@data[i] <- log_ens@data[i]+10*log10(r$layer)
        ens_nc@data[i] <- 10^(log_ens@data[i]/10)
      }
      
    } else { 
      
      # make a spatial point object out of the observations
      names(err) <- c("err", "X", "Y")
      coordinates(err) <- ~X+Y
      err@proj4string@projargs <- myproj
      err <- crop(err, e)
      
      if (is.null(err)==T){ # in case the observation variable is not empty but it's all NA the simulation is still unconditional
        
        #define prediction grid
        datagrid <- r
        datagrid$layer <- datagrid$layer*0
        datagrid@proj4string@projargs <- myproj
        
        # unconditional simulation
        g.obj <- gstat(formula=z~1, locations = ~x+y, dummy=T, beta=mu, model=variog, nmax=24)
        sim <- predict(g.obj, datagrid, nsim = 100)
        
        # Generate the ensembles
        log_ens <- sim
        ens_nc <- sim
        for (i in 1:100){
          log_ens@data[i] <- log_ens@data[i]+10*log10(r$layer)
          ens_nc@data[i] <- 10^(log_ens@data[i]/10)
        }
        
      } else { # case in which we have observations for a conditional simulation
        
        #define prediction grid
        datagrid <- r
        datagrid$layer <- datagrid$layer*0
        datagrid@proj4string@projargs <- myproj
        
        # conditional simulation according to Delhomme 1979
        
        # 1) non conditional simulation
        g.obj <- gstat(formula=z~1, locations = ~x+y, dummy=T, beta=mu, model=variog, nmax=24)
        sim_un <- predict(g.obj, datagrid, nsim = 100)
        
        # 2) kriged errors
        Krigerr <- krige(err~1, err, newdata=datagrid, model=variog)
        
        # 3) kriged simulations at error locations and calculate conditional simulation
        extr_un <- over(err, sim_un)
        sim <- sim_un
        for (i in 1:100){
          cond <- err
          cond$err <- unlist(extr_un[i])
          kcond <- krige(err~1, cond, newdata=datagrid, model=variog)
          sim@data[i] <- Krigerr$var1.pred + sim_un@data[i] - kcond$var1.pred
        }
        
        # generate the ensembles
        log_ens <- sim
        ens_nc <- sim
        for (i in 1:100){
          log_ens@data[i] <- log_ens@data[i]+10*log10(r$layer)
          ens_nc@data[i] <- 10^(log_ens@data[i]/10)
        }
      }
    }
    
    
    # Ensemble correction
     
    g <- as.numeric(G0[t,])
    g <- data.frame(g,coord)
    g <- g[complete.cases(g),]
    names(g) <- c("rg", "X", "Y")
    coordinates(g) <- ~X+Y
    g@proj4string@projargs <- myproj
    g <- crop(g, e)
    
    sigma_g <- sd(g$rg, na.rm=T)
    mu_g <- mean(g$rg, na.rm=T)
    
    if (sigma_g!=0 & mu_g!=0){
      ens_extr <- over(g, ens_nc)
      ens <- ens_nc
      sigma_en <- rep(0,100)
      mu_en <- rep(0,100)
      for (i in 1:100){
        sigma_en <- sd(unlist(ens_extr[i]), na.rm=T)
        mu_en <- mean(unlist(ens_extr[i]), na.rm=T)
      }
      sigma_en <- mean(sigma_en, na.rm=T)
      mu_en <- mean(mu_en, na.rm=T)
      for (i in 1:100){
        ens@data[i] <- (sigma_g/sigma_en)*(ens_nc@data[i]-mu_en)+mu_g
        ens@data[(r$layer<0.1),i] <- r$layer[(r$layer<0.1)]
      }
    } else {
      ens <- ens_nc
    }
    ens@data[ens@data<0] <- 0
    if(dates[t,1]==2008 & dates[t,2]==1 & dates[t,3]==1 & dates[t,4]==10){
      v <- c(nugget_spec,sill_spec,range_spec)
      save(r, ens, ens_nc, g, v, file=paste0(results_path, "ensemble-t", substr(radarname,6,17)))
    }
    
    # extract on validation RG
    
    g_val <- as.numeric(G_val[t,])
    g_val <- data.frame(g_val,coord_val)
    g_val <- g_val[complete.cases(g_val),]
    names(g_val) <- c("rg", "X", "Y")
    coordinates(g_val) <- ~X+Y
    g_val@proj4string@projargs <- myproj
    
    for (i in 1:100){
      ens_val[,i] <- unlist(over(g_val, ens[i]))
    }
    
    # Extract on 3 basins
    
    # Extract ensembles on three basins
    ens_lune <- over(Lune, ens, FUN="mean")
    ens_rawthey <- over(Rawthey, ens, FUN="mean")
    ens_ribble <- over(Ribble, ens, FUN="mean")
    
    # Save into hydro variable
    ens_hydro[1,] <- as.numeric(ens_lune)
    ens_hydro[2,] <- as.numeric(ens_rawthey)
    ens_hydro[3,] <- as.numeric(ens_ribble)
    
    save(nugget_spec, sill_spec, range_spec, n_hours, file=paste0(variog_path, "Variogram_parameters_", toString(t)))
    save(ens_val, file=paste0(val_path, "ens_val_", toString(t)))
    save(ens_hydro, file=paste0(hydro_path, "ens_hydro_",toString(t)))
    
  }
}



