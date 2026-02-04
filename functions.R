

# Libraries ---------------------------------------------------------------

library(ncdf4) 
library(abind)
library(lubridate)
library(kSamples)
library(spatstat.univar)

# Read data function ------------------------------------------------------

read.data <- function(dir, index, ssp, model, ens.members=NULL){
  
  print(paste0("Reading in data for ", model$name, ", ", index, " and SSP", ssp))

  # create string of ensemble members to read
  if (length(ens.members) != 0 & length(ens.members)!= model$members) stop("List of ensemble members is too long or too short")
  if (length(ens.members) != 0){
    ens.members <- paste0("(", paste(ens.members, collapse="|"), ")")
  }
  
  ## list files for historical period
  f.hist <- list.files(dir, pattern=paste0("^", index, ".*yr_", model$name, "_historical_", ens.members), full.names=T)
  if (length(f.hist)!=model$members) stop("Historical period: incorrect number of ensemble members")
  
  ## list files for future period
  f.ssp <- list.files(dir, pattern=paste0("^", index, ".*yr_", model$name, "_ssp", ssp, "_", ens.members), full.names=T)
  if (length(f.ssp)!=model$members) stop("Future period: incorrect number of ensemble members")

  ## loop over files (ensemble members) and combine into an array
  data.all <- NULL
  for (f in seq_along(f.hist)){
    
    ## read and check historical data
    foo.hist <- nc_open(f.hist[f])  
      var.hist <- ncvar_get(foo.hist, paste0(index, "ETCCDI"))
      times.hist <- as.POSIXct(foo.hist$dim$time$vals*3600*24, origin="1850-01-01", tz="Etc/UTC")
    nc_close(foo.hist)
    if (length(times.hist) < 165) stop("Historical period too short")
    
    ## read and check future data
    foo.ssp <- nc_open(f.ssp[f])
      var.ssp <- ncvar_get(foo.ssp, paste0(index, "ETCCDI"))
      times.ssp <- as.POSIXct(foo.ssp$dim$time$vals*3600*24, origin="1850-01-01", tz="Etc/UTC")
    nc_close(foo.ssp)
    if (length(times.ssp) > 86){
      ind <- (year(times.ssp) %in% 2015:2100)
      times.ssp <- times.ssp[ind]
      var.ssp <- var.ssp[,,ind]
    }
    if (length(times.ssp) < 86) stop("Future period too short")
    
    ## Combine data
    var.foo <- abind(var.hist, var.ssp, along=3)
    data.all <- abind(data.all, var.foo, rev.along=4)
    dimnames(data.all)[[4]] <- year(c(times.hist, times.ssp))

  }
  return(data.all)
}
