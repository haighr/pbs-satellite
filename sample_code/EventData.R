library(PBSmapping)
## ncdf4 needed for reading netCDF data
library(ncdf4)

nceventdata_get <- function(filename, week=1, dataname="DATA",proj="DATA"){
    nc <- nc_open(filename, readunlim=FALSE)

    ## obtain the measurements and coordinates from the NetCDF file
    dataVariable <- nc$var[[1]]
    data <- ncvar_get(nc, dataVariable)
    lon <- ncvar_get(nc, "lon")
    lat <- ncvar_get(nc, "lat")

    nc_close(nc)
    
    ## vector length for one week of data
    size <- length(lon) * length(lat)

    ## build the data frame
    df <- data.frame(EID=1:size)
    df$X <- rep(lon, times=length(lat))
    df$Y <- rep(lat, each=length(lon))

    offset <- (week - 1) * size
    df[[dataname]] <- data[offset:(offset + size)]

    df <- as.EventData(df, projection=proj)

    return(df)
}

