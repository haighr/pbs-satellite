library(ncdf4)


read.ncdfData <- function (filename,  dataVariable=1, 
                           convertMissingValues=FALSE, dataType=NULL,
                           dataUnits=NULL,
                           xlim=NULL, ylim=NULL, tlim=NULL,
                           x=NULL,y=NULL, time=NULL,
                           Ux=NULL, Uy=NULL, Utime=NULL) {
    
    nc <- nc_open(filename, readunlim=FALSE)
    
    ## get the data object specified by the user
    data <- nc$var[[dataVariable]]

    ## if unspecified by the user (NULL), set sensible defaults
    if (is.null(dataType)) {
        dataType <- data$longname
    }
    
    if (is.null(dataUnits)) {
        dataUnits <- data$units
    }
    
    if (is.null(time)) {
        time <- "time"
    }
    
    if (is.null(x)) {
        x <- "lon"
    }
    
    if (is.null(y)) {
        y <- "lat"
    }
    
    ## anytime the user did not specify a dimension's unit, attempt to
    ## extract it automatically
    for (n in 1:length(data$dim)) {
        name <- data$dim[[n]]$name
        if (is.null(Ux) && name==x) {
            Ux <- data$dim[[n]]$units
        }
        else if (is.null(Uy) && name==y) {
            Uy <- data$dim[[n]]$units
        }
        else if (is.null(Utime) && name==time) {
            Utime <- data$dim[[n]]$units
        }
    }
    
    ## if units not found, we also could not find the dimension
    if (is.null(Ux)) stop("x dimension not found; specify using 'x'")
    if (is.null(Uy)) stop("y dimension not found; specify using 'y'")
    if (is.null(Utime)) stop("time dimension not found; specify using 'time'")
    
    ## store x and y vectors containing the spatial representation of data
    ## points
    xvect <- as.vector(ncvar_get(nc,x))
    yvect <- as.vector(ncvar_get(nc,y))
    
    ## store length of x and y vectors
    xlen <- length(xvect)
    ylen <- length(yvect)
    
    ## split up time units into a vector
    ## e.g., time units might be "days since 1900-01-01 00:00:00"
    timeV <- strsplit(Utime, " ")[[1]]
    
    ## get the actual time stamps from nc file
    names <- ncvar_get(nc, time)
    
    ## convert time stamps to strings in the format
    ## "YYYY-MM-DD HH:MM:SS"
    names <- .getSliceNames(names, unit=timeV[1],
                            origin=paste(timeV[3], timeV[4]))
    
    ## user has specified a tlim argument, need to limit the amount of slices
    if (!is.null(tlim)) {

        ## TO DO: force the user to provide the time as a string -- then we can
        ## ensure the conversion happens with tz="GMT" argument
        if (length(tlim) != 2) {
            nc_close(nc)
            stop("tlim must be a vector of two strings, each formatted as \"YYYY-MM-DD HH:MM:SS\"")
        }

        ## convert date strings into POSIXlt objects for date comparison
        names <- as.POSIXlt(names, tz="GMT", origin=paste(timeV[3], timeV[4]))

        ## keep only the dates that fall in tlim
        datesInTLIM <- names >= as.POSIXlt(tlim[1], tz="GMT") &
            names <= as.POSIXlt(tlim[2], tz="GMT")

        ## get locations of the data sheets in the data object
        rangeT <- which(datesInTLIM)
        
        ## use locations of first and last value in rangeT to create of offset into the data object
        ## this is done to prevent unnecessary processing of data that is not in tlim
        dataset <- ncvar_get(nc, data,
                             start=c(1, 1, min(rangeT)),
                             count=c(xlen, ylen, max(rangeT)-min(rangeT)+1))
        
        ## convert POSIXlt objects back to strings for name assignment
        names <- as.character(names[datesInTLIM])
        
        ## no tlim specified, read in the full contents of data object  
    } else { 
        dataset <- ncvar_get(nc, data)
    }

    ## clip the data to keep only data within xlim/ylim
    if (!is.null(xlim)) {
        xvect <- xvect[xvect >= xlim[1] & xvect <= xlim[2]]
        dataset <- dataset[xvect,,]
        xlen <- length(xvect)
    }
    
    if (!is.null(ylim)) {
        yvect <- yvect[yvect >= ylim[1] & yvect <= ylim[2]]
        dataset <- dataset[,yvect,]
        ylen <- length(yvect)
    }

    ## TO DO: how will you change this? replacements are automatically
    ## made...
    if (convertMissingValues == TRUE) {
        missval <- data$missval
        missingLocations <- which(dataset == missval)
        dataset[missingLocations] <- NA
    }
    
    ## calculate size of each "data sheet" in the nc file
    size <- (xlen * ylen)
    
    ## calculate number of data sheets in the dataset
    length <- length(dataset) / size
    
    ## convert dataset into a list of vectors
    ncdfData <- split(dataset, rep(1:length, each=size))
    
    ## put dimension on each list vector to convert to a list of matrices;
    ncdfData <- lapply(ncdfData, ".createDim", xlen, ylen)
    
    ## apply name changes to the list 
    names(ncdfData) <- names
    
    ## add attributes to ncdfData object
    attr(ncdfData,"dataType") <- dataType
    attr(ncdfData,"dataUnits") <- dataUnits
    attr(ncdfData,"x") <- xvect
    attr(ncdfData,"y") <- yvect
    ## TO DO: make the value of this attribute consistent with
    ## with whatever the data set is presently using to represent missing data
    attr(ncdfData,"missval") <- data$missval
    
    nc_close(nc)
    
    return(ncdfData)
}

## function applies a dimension to a list of vectors
.createDim <- function(vectors, xlen, ylen){
    ## at this point, one would index the first dimension with
    ## an x value and the second with a y value -- in the
    ## order you would expect for indexing, i.e., (x, y)
    dim(vectors) <- c(xlen, ylen)
    
    return(vectors)
}

## function converts a vector of time stamps to YYYY-MM-DD HH:MM:SS format 
.getSliceNames <- function(times, unit, origin) {
    if(unit=="days"){
        times <- as.character(as.POSIXlt(times * (24 * 60 * 60), tz="GMT", origin=origin))
    } else if(unit=="seconds"){
        times <- as.character(as.POSIXlt(times, tz="GMT", origin=origin))
    } else if(unit=="hours"){
        times <- as.character(as.POSIXlt(times  * (60 * 60), tz="GMT", origin=origin))
    } else {
        warning("Was not able to properly label slice names")
    }
    return (times)
}


time <- c("1999-01-01","2000-01-01")
time <- as.POSIXlt(time, tz="GMT",format="%Y-%m-%d")

l <- read.ncdfData("~/Desktop/nc/sst.wkmean.1990-present.nc",tlim=time,xlim=c(50,100),ylim=c(50,100))

l <- read.ncdfData("~/Desktop/nc/chlonasa19972002.1-deg.nc")

l <- read.ncdfData("/home/nicholas/Desktop/nc/19811231-NCDC-L4LRblend-GLOB-v01-fv02_0-AVHRR_OI.nc")

l <- read.ncdfData("~/Desktop/nc/CRUTEM.4.3.0.0.anomalies.nc",x="longitude",y="latitude")

l <- read.ncdfData("~/Desktop/nc/LASoutput.nc",x="longitude",y="latitude")

l <- read.ncdfData("~/Desktop/nc/LASoutput.nc",x="LONGITUDE",y="LATITUDE",dataVariable=2)

l <- read.ncdfData("~/Desktop/nc/LASoutput.nc",x="LONGITUDE",y="LATITUDE",time="TIME",dataVariable=2)

names(l)


