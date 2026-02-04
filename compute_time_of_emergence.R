

# Functions and namelists ---------------------------------------------------------------

source("functions.R")

index.list <- c("altcdd", "csdi", "altcwd", "dtr", "fd", "gsl", "id", "prcptot", "r10mm", "r1mm", "r20mm", "r95p", "r99p", "rx1day", "rx5day", "sdii", "su", "tn10p", "tn90p", "tnn", "tnx", "tr", "tx10p", "tx90p", "txn", "txx", "wsdi")

ssp.list <- c(126, 245, 370, 585)

model.list <- data.frame(
  name=c("ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3", "EC-Earth3-Veg", "FGOALS-g3", "GFDL-ESM4", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "MIROC6", "MIROC-ES2L", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM", "NorESM2-MM", "UKESM1-0-LL"), 
  members=c(1, 40, 1, 50, 1, 1, 7, 1, 1, 1, 1, 1, 1, 3, 1, 1, 30, 1, 1, 1, 1))



# Example: one model, index and scenario ----------------------

## directory containing files for annual ETCCDI indices
## files are of the format "tnnETCCDI_yr_ACCESS-ESM1-5_historical_r10i1p1f1_1850-2014.nc"

dir.ind <- "/div/amoc/ClimateIndicatorDB/Global/raw/hist/tempCMIP6/Climdex_base1981-2010/"

i <- 20  # tnn
s <- 4   # SSP5-8.5
m <- 14   # MIROC6
ens.members <- c("r1i1p1f1", "r2i1p1f1", "r3i1p1f1")


# Read data ---------------------------------------------------------------

# longitude is ordered from -180 to +180, latitude from -90 to + 90
data.all <- read.data(dir.ind, index.list[i], ssp.list[s], model.list[m,], ens.members=ens.members)

## read in land-sea mask
file.lsm <- list.files(dir.ind, pattern=paste0("sftlf.*", model.list$name))
foo <- nc_open(file.lsm)
  lsm <- ncvar_get(foo, "sftlf")
  lon.val <- foo$dim$lon$vals
  lat.val <- foo$dim$lat$vals
nc_close(foo)


# Run KS test and store the p values --------------------------------------

years <- 1850:2100
periods <- 1901:2081

data.ref <- data.all[,,,as.character(1850:1900),drop=FALSE] # reference period
ks.values <- array(NA, dim=c(model.list$members[m], length(periods), length(lon.val), length(lat.val)))

for (p in seq_along(periods)){ # for every 20-year period from 1901-2081
  
  print(paste0(periods[p], "-", periods[p]+19))
  
  years.p <- seq(periods[p], periods[p]+19)
  data.period <- data.all[,,,as.character(years.p),drop=FALSE]
  
  for (lon in seq_along(lon.val)){ # for every grid point
    for (lat in seq_along(lat.val)){
      
      if (lsm[lon,lat]<=50) next # skip non-land grid points
      
      for (e in 1:model.list$members[m]){ # for every ensemble member
      
        data1 <- data.ref[e,lon,lat,,drop=FALSE]
        data2 <- data.period[e,lon,lat,,drop=FALSE]
        
        ks.values[e,p,lon,lat] <- ks.test(data1,data2)$p.value
      } 
    }
  }
}



# Find time of emergence based on grid-wise p-values ----------------------

emergence.lonlat <- array(NA, dim=c(model.list$members[m], length(lon.val), length(lat.val)))
  
ind.TF <- (ks.values<0.05) # p-values below significance level?

for (lon in seq_along(lon.val)){
  for (lat in seq_along(lat.val)){
    for (e in 1:model.list$members[m]){
      
      TF.ts <- ind.TF[e,,lon,lat] # time series for grid point
      if (all(is.na(TF.ts))) next # skip non-land grid points
            
      for (p in seq_along(periods)){
        
        if (is.na(TF.ts[p])) next # skip period if value is NA
        
        if (p > 170){ # if period has reached 2070, mark grid point as no emergence (here a large value)
          emergence.lonlat[e,lon,lat] <- 2100
        break
        }
        
        if (!TF.ts[p]) next # if there is no significance in this time period, move to the next
        
        TF.tail <- tail(TF.ts, -p) # remaining part of the time series
        if (sum(TF.tail, na.rm=T) > 0.95 * length(TF.tail)){ # if 95% of the remaining periods are significant, record period as time of emergence
          emergence.lonlat[e,lon,lat] <- periods[p]
          break
        }
      }
    }
  }
}


# Aggregate over regions (median) --------------------------------------

## load file with region info
file.regions <- paste0(dir.ind, "regions", model.list$name[m], "_AR6.Rdata")
load(file.regions) # longitude is ordered from -180 to +180, latitude from -90 to + 90

## define land regions
seas <- sort(unique(names[types=="Ocean"]))
region.names <- labels[!labels$name %in% seas,]
region.names <- region.names[order(region.names$name),]

emergence.reg.median <- array(NA, dim=c(model.list$members[m], nrow(region.names)))
  
for (e in 1:model.list$members[m]){
  for (r in 1:nrow(region.names)){
    ind.r <- which(names==region.names$name[r], arr.ind=T) # which grid points belong to region
    data.r <- mapply(function(k,l) emergence.lonlat[e,k,l], ind.r[,1], ind.r[,2]) # subset these grid points
    emergence.reg.median[e,r] <- median(data.r, na.rm=T)
  }
}


# Aggregate over models (weighted median) ---------------------------------

## run code above for all models first, store data in emergence.reg.median in a list all.models.region

## load file with weights info
tas <- TRUE # TRUE for temperature indices, FALSE for precipitation indices
file.weights <- paste0(dir.ind, "weights_", ifelse(tas, "tas", "pr"), ".Rdata")
load(file.weights)

if (tas) { # remove weight for NorESM2-LM for temperature and redistribute
  weights["NorESM2-LM"] <- 0
  weights <- weights/sum(weights)
}

## ToE
emergence.median <- drop(abind(all.models.region, along=1)) # collapse list to an array of dim c(length(weights), nrow(region.names))
model.w.median <- apply(emergence.median, 2, weighted.median, w=weights) # weighted median across regions, values > 2070 correspond to no emergence

## model agreement
model.median <- t(sapply(all.models.region, function(x){apply(x, 2, median, na.rm=T)})) # median across individual ensemble members for each model
model.agree <- apply(model.median<=2070, 2, sum, na.rm=T) # sum of models with emergence before 2070









