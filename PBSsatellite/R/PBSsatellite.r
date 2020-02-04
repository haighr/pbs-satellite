##==============================================================================
## PBSsatellite
## ------------
##  addVstrip.................Add a vertical colour strip as a legend
##  assessMissingData.........Calculate percent missing data from each slice
##  clipRegion................Clip a polygon to a smaller region
##  convert.ncdfData..........Convert an ASCII file to a netCDF (ncdf4) file
##  create.ncdfData...........Create a NetCDF (ncdf4) file from data
##  download.ncdfData.........Download NetCDF satellite data files
##  extractSlices.............Extract slices from an ncdfData object
##  extractTimeSeries.........Extract time series from an ncdfData object
##  listToDF..................Convert a list to a data.frame
##  merge.ncdfData............Merge a suite of netCDF binary files into one ncdfData
##  plot.ncdfData.............Plot an ncdfData slice
##  print.ncdfData............Print a summary of the ncdfData object
##  read.ncdfData.............Read netCDF binary file and create ncdfData object in R
##  removeAnomalousValues.....Remove anomalous values from ncdfData
##  scaleRegion...............Scale ncdfData slices to a new resolution based on a scale factor
##  to.EventData..............Convert an ncdfData object to an EventData object
##  write.ncdfData............Write slices of ncdfData to CSV files
##
##-----Supplementary hidden functions-----
##  .createDim................Creates matrix from list of vectors then puts new matrix into a data layer
##  .createScaledDownSlices...Helper function to scale down slices
##  .createScaledUpSlices.....Helper function to scale up slices
##  .createScaleSplitVector...Creates vector of indices used to extract data from ncdfData slice
##  .findRC...................Return no. (rows, columns) for multi-panel figures given no. figures to fit on one page
##  .getClippedMatrix.........Mask of T|F values of minimum matrix size to fit full size of polygon
##  .getFunctionResults.......Applies function(s) to list of sliceValues
##  .getNumberSequenceInfo....Returns a list containing numberSequence vector used for clipping vectors
##  .getPolysVectorList.......Returns a list detailing which polygons vertices fall inside a polygon region
##  .getSlice.................Extracts a slice from an ncdfData object
##  .getSliceNames............Converts vector of time stamps to YYYY-MM-DD HH:MM:SS format
##  .ncdfDataTlimClip.........Removes data from ncdfData that does not fall in the specified tlim period
##  .ncdfDataXClip............Removes data from ncdfData that does not fall in the specified xlim region
##  .ncdfDataYClip............Removes data from ncdfData that does not fall in the specified ylim region
##  .pInv.....................Find nearest position in vector choice using a target point
##===============================================================================


## addVstrip----------------------------2018-11-08
## Add a vertical colour strip as a legend.
## ---------------------------------------------RH
addVstrip = function (x, y, zlim, col, title, xwidth=0.02, yheight=0.5, ...) 
{
	#if (dev.cur()>1) { oldpar=par(no.readonly=TRUE); on.exit(par(oldpar)) }
	uxy <- par()$usr
	x1 <- uxy[1];  x2 <- uxy[2]
	y1 <- uxy[3];  y2 <- uxy[4]
	x0 <- x1 + x * (x2 - x1)
	y0 <- y1 + y * (y2 - y1)
	xw0 = xwidth * (x2-x1)
	yh0 = yheight * (y2-y1)
	if (par()$xlog){
		x0 <- 10^x0; xw0 = 10^xw0 }
	if (par()$ylog){
		y0 <- 10^y0; yh0 = 10^yh0 }
	xval = x0 + c(0,xw0)
	yval = seq(y0-yh0,y0,len=length(col)+1)
	ncol = length(col)
	xpol = c(x0,x0,xval[2],xval[2],NA)
	xpol = rep(xpol,ncol)
	ypol = numeric()
	for (i in 1:ncol)
		ypol = c(ypol, c(yval[i],rep(yval[i+1],2),yval[i],NA))
	polygon(xpol,ypol,col=col, border=FALSE)

	## Strip outline
	xbox = xval[c(1,1,2,2,1)]
	ybox = c(y0-yh0,y0)[c(1,2,2,1,1)]
	lines(xbox, ybox, lwd=0.5, col="grey30")

	## Interpolate break points for labelling
	ybrks = sort(unique(ypol))
	zbrks = seq(zlim[1], zlim[2], length.out=length(ybrks))
	zuse  = pretty(zlim, n=10)
	zuse  = zuse[zuse>min(zbrks) & zuse<max(zbrks)]
	yuse  = approx(zbrks,ybrks,zuse)$y
	nbrks = length(yuse)
	xuse  = rep(xval[2],nbrks)
	xlin  = rep(c(xval,NA), nbrks)
	ylin  = as.vector(sapply(yuse,function(x){c(x,x,NA)}))
	lines(xlin, ylin, lwd=0.5, col="grey30")
	text(xuse, yuse, zuse, pos=4, offset=0.25, cex=0.8)
	if (!missing(title))
		text(mean(xval), max(ybox), title, pos=3, offset=0.5)
	invisible(list(xbox=xbox, ybox=ybox, xuse=xuse, yuse=yuse, zuse=zuse))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~addVstrip


## assessMissingData--------------------2016-03-23
##  Calculate percent missing data from each slice
## ------------------------------------------NL/RH
assessMissingData <- function(ncdfData, tlim=NULL, xlim=NULL, 
   ylim=NULL, polygons=NULL, include.lowest=TRUE) {
	## validate ncdfData is a ncdfData object
	if (class(ncdfData) != "ncdfData") {
		stop("ncdfData must be of class ncdfData")
	}
	## if user provide a tlim argument
	if(!is.null(tlim)){
		ncdfData <- extractSlices(ncdfData, tlim)
	}
	## provide clipping if necessary 
	if(!is.null(xlim) || !is.null(ylim) || !is.null(polygon)){
	ncdfData <- clipRegion(ncdfData, xlim=xlim, ylim=ylim,
		polygons=polygons,
		include.lowest=include.lowest)
	}
	## creates a list of missing data percentages from each slices data layer
	missingVectorsList <- lapply(ncdfData, function(slice)
	{
		## ncdfData is read in to be a list of lists, where each list is called "data" and contains a single spatial matrix.
		## This is differnt from a list of matrices.
		#nanCount <- sum(is.nan(slice$data))  ## does not compute on a list
		#return(((length(slice$data[is.na(slice$data)]))-nanCount/length(slice$data))*100) ## makes no sense whatsoever (RH)
		return( 100. * sum(is.na(slice$data)) / length(slice$data) )
	})
	## unlist to get a vector of missing percentages
	return(unlist(missingVectorsList))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~assessMissingData


## clipRegion---------------------------2018-10-31
## Clip a polygon to a smaller region
## ---------------------------------------------NL
clipRegion <- function(ncdfData, xlim=NULL, ylim=NULL,
   polygons=NULL, include.lowest=TRUE)
{
	## save original ncdfData attributes in case clipping by polygon yield one X and/or Y
	ncdfDataAttributes.original <- attributes(ncdfData)

	## validate ncdfData is a ncdfData object
	if (class(ncdfData) != "ncdfData") {
		stop("ncdfData must be of class ncdfData")
	}
	## if user has provided an xlim agrument 
	if(!is.null(xlim)) {
		ncdfData <- .ncdfDataXClip(ncdfData, xlim, include.lowest)
	}
	## if user has provided a ylim argument
	
	if(!is.null(ylim)) {
		ncdfData <- .ncdfDataYClip(ncdfData, ylim, include.lowest)
	}
	## if the user has specified polygons
	if (!is.null(polygons)) {
		## retrieve x and y vector attributes from ncdfData
		xvect <- attr(ncdfData, "x")
		yvect <- attr(ncdfData, "y")

		## convert first slice to EventData
		## note: only need to clip one slice and then apply clipping result
		## to all existing slices
		ed <- to.EventData(ncdfData, 1)

		## get the contents of the matrix mask with the xvect and yvect sized
		## to insure matrix is square while fitting the polygons in its entirety
		## NaN values are used to represent pixels of clipped data
		clipResult <- .getClippedMatrix(ed, xvect, yvect, polygons)

		## when clipping to polygons the coordinate span of a vector also changes
		## in most cases
		newXVect <- xvect[clipResult$xIDx]
		newYVect <- yvect[clipResult$yIDx]

		## location in xvect where new xlims are  
		minX <- min(clipResult$xIDx)
		maxX <- max(clipResult$xIDx)
		
		## location in yvect where new ylims are 
		minY <- min(clipResult$yIDx)
		maxY <- max(clipResult$yIDx)
		
		## assign new x and y coordinate vectors to fit polygons in range
		attr(ncdfData, "y") <- newYVect
		attr(ncdfData, "x") <- newXVect
		if (length(newXVect)==1)
			attr(ncdfData, "xBy") <- abs(unique(round(diff(ncdfDataAttributes.original$x), 3)))
		if (length(newYVect)==1)
			attr(ncdfData, "yBy") <- abs(unique(round(diff(ncdfDataAttributes.original$y), 3)))

		## make backup of ncdfData attributes
		ncdfDataAttributes <- attributes(ncdfData)

		## applies slice mask to all layers in a slice
		## regions that do not fall in the polygon/slice mask will be replaced with NaN
		## a slices x and y vectors will be cropped down to have the minimal space region
		## in order to fully store all polygons in a matrix

		ncdfData <- lapply(ncdfData, function(slice, xlim, ylim, sliceMask)
		{
			for(layerName in names(slice)){
				## gets the minimal region of data to fit all desired polygons
				slice[[layerName]] <- slice[[layerName]][xlim[1]:xlim[2], ylim[1]:ylim[2], drop=FALSE] 
				## converts everything outside of sliceMask/polygon regions to NaN
				slice[[layerName]][which(!sliceMask)] <- NaN
			}
			return(slice)
		}, xlim=c(minX,maxX), ylim=c(minY,maxY), sliceMask=clipResult$sliceMask)

		## reassign ncdfData attributes
		attributes(ncdfData) <- ncdfDataAttributes
	}
	return(ncdfData)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~clipRegion


## convert.ncdfData---------------------2016-03-22
## Convert an ASCII file to a netCDF (ncdf4) file.
## ---------------------------------------------RH
convert.ncdfData = function(filename, zfld, nc.filename="converted.nc",
   summary.func=sum, offset=c(0,0), mv=-99, 
   dataType="Chlorophyll", dataUnits="mg/m3")
{
	errmess = as.character()
	textfile = fread(filename)
	class(textfile) = "data.frame"

	flds = names(textfile)
	xfld = flds[grep("^X$|[Ll]on",flds)][1]
	yfld = flds[grep("^Y$|[Ll]at",flds)][1]
	tfld = flds[grep("^T$|[Dd]ate",flds)][1]
	if(any(is.na(c(xfld,yfld,tfld))))
		stop( switch(c(1:3)[is.na(c(xfld,yfld,tfld))][1],
		"xfld missing -- need `X` or some variation on `longitude`",
		"yfld missing -- need `Y` or some variation on `latitude`",
		"tfld missing -- need `T` or some variation on `date`") )
	if (!any(zfld %in% flds))
		stop("zfld must specify at least one field name (e.g.,`Chl`) in the data file")

	textfile[[xfld]] = textfile[[xfld]] + offset[1]
	textfile[[yfld]] = textfile[[yfld]] + offset[2]

	xvals = sort(unique(textfile[[xfld]])); nx = length(xvals)
	yvals = sort(unique(textfile[[yfld]])); ny = length(yvals)
	Tvals = sort(unique(textfile[[tfld]])); nt = length(Tvals)

	tvals = as.numeric(as.Date(Tvals) - as.Date(Tvals[1]))
	tmess = paste0("days since ",Tvals[1])

	zlist = list()
	for (z in zfld) {
		if (!(z%in%flds)) next 
		zmat  = array(NA, dim=c(nx,ny,nt), dimnames=list(x=xvals,y=yvals,t=Tvals))
		zfile = split(textfile,textfile[[tfld]])
		for (tt in names(zfile)) {
			tfile = zfile[[tt]]
			tfile[[z]][is.element(tfile[[z]],mv)] = NA
			if (all(is.na(tfile[[z]]))) next
			expr = paste0("afile = aggregate(", z ," ~ ", xfld, " + ", yfld, ", data=tfile, FUN=",substitute(summary.func),"); ")
			expr = paste0(expr, "amat = xtabs(",z,"~",xfld,"+",yfld,", data=afile)")
			eval(parse(text=expr))
			amat[amat==0] = NA
			zmat[dimnames(amat)[[1]], dimnames(amat)[[2]], tt] = amat
		}
		zlist[[z]] = zmat
	}
	create.ncdfData(filename=nc.filename, xvals=xvals, yvals=yvals, tvals=tvals, tmess=tmess, 
		dataType=dataType, dataUnits=dataUnits, zlist=zlist)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~convert.ncdfData


## create.ncdfData----------------------2016-03-21
## Create an netCDF (ncdf4) file from data.
## See posting by by user3710546  (Mar 11 '15 at 10:48) at:
## http://stackoverflow.com/questions/28949971/writing-data-to-a-netcdf-file-with-r
## ---------------------------------------------RH
create.ncdfData = function(filename, xvals, yvals, tvals, 
   tmess="days since 1900-01-01", zlist, mv=-99, 
   dataType="Chlorophyll", dataUnits="mg/m3", longname=filename)
{
	nx <- length(xvals)
	ny <- length(yvals)
	nt <- length(tvals)

	lon  <- ncdim_def("lon", "degrees_east", xvals)
	lat  <- ncdim_def("lat", "degrees_north", yvals)
	time <- ncdim_def("time",tmess, tvals, unlim=TRUE)

	dType = dataType; dUnit = dataUnits
	if (length(dataType)<length(zlist)) {
		dType = rep(dataType,length(zlist))[1:length(dataType)]
		dType = paste0("z",1:length(zlist),".",dType)
		dUnit = rep(dataUnits,length(zlist))[1:length(dataUnits)]
		dUnit = paste0("z",1:length(zlist),".",dUnit)
	}
	## The source file may have more than one column of z-data
	defList = list()
	for (z in 1:length(zlist)) {
		zz = names(zlist)[z]
		defList[[paste0(zz,"_def")]] = ncvar_def(dType, dUnit, list(lon, lat, time), longname=longname, mv)
	}
	ncnew <- nc_create(filename, vars=defList)
	for (z in 1:length(zlist)) {
		ncvar_put(ncnew,defList[[z]],zlist[[z]])
	}
	nc_close(ncnew)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~create.ncdfData


## download.ncdfData--------------------2018-11-05
## Download NetCDF satellite data files.
## ---------------------------------------------RH
download.ncdfData <- function(url, save.dir="./rawdata",
   years, months, days, res, pname, period, sensor, curlbin,
   need.seq=FALSE, overwrite=FALSE, verbose=FALSE)
{
	## Get rid of trailing delimiters (if they are specified)
	url      = sub("/$","",url)
	save.dir = sub("/$","",save.dir)
	if (!missing(curlbin))
		curlbin = sub("/$","",curlbin)

	.flush.cat = function (...) {
		cat(...); flush.console(); invisible() }

	if (need.seq) {
		if(!missing(years) && length(years) != length(min(years):max(years)))
			stop("years vector is not a sequence ascending by 1")
		if(!missing(months) && length(months) != length(min(months):max(months)))
			stop("months vector is not a sequence ascending by 1")
		if(!missing(days) && length(days) != length(min(days):max(days)))
			stop("days vector is not a sequence ascending by 1")
	}
	if(!dir.exists(save.dir)) {
			dir.create(save.dir, recursive = TRUE)
	}
	surl = strsplit(url, split="/")[[1]]

	## -----------------
	## ERSST.v4 SST data
	## -----------------
	if (grepl("ersst",url))
	{
		mess = c(
		"\nAdapted from Michael Malick's revised 'sst_download' function.",
		"Download monthly ERSST.v4 SST data in netCDF format",
		"https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v4/netcdf",
		"-----",
		"Filename convention is:",
		"   ersst.VERSION.yyyymm.nc",
		"   yyyy=four digit year",
		"   mm=two digit month",
		"-----",
		"years  = numeric vector of years to download data for, e.g., 1950:2016",
		"         should be a numeric sequence, ascending by 1",
		"months = numeric vector of months to download data, e.g., 1:12",
		"         should be a numeric sequence, ascending by 1",
		"save.dir = directory to save sst data files in, e.g., './data/'"
		)
		if (verbose)
			.flush.cat(paste0(mess,collapse="\n"),"\n\n")
		if(missing(years) || missing(months))
			stop("Supply 'years' and 'month' for this URL.")
		vers = surl[grep("v[0-9]", surl)]
		cnt  = 0
		for(i in years) {
			for(j in months) {
				if(j < 10)
					j  = paste0("0", j)
				fname = paste0("ersst", ifelse(vers=="v3b","",paste0(".",vers)), ".", i, j, ".nc")
				web   = paste0(url, "/", fname)
				fpath = paste0(save.dir, "/", fname)
				if(!file.exists(fpath) || overwrite) {
					download.file(web, destfile=fpath, method="auto", mode="wb")
					cnt = cnt + 1
				}
			}
		}
		## Check if all data files were downloaded
		n.possible <- length(years) * length(months)
	}
	## ----------------
	## Ocean Color Data
	## ----------------
	if (grepl("oceandata",url)) 
	{
		mess = c(
		"\nThe Ocean Color (OC) products all include combinations of the following derived parameters (pname in parentheses):",
		"  nLw    -- normalized water-leaving radiance (???);",
		"  Rrs    -- remote sensing reflectance (RRS_Rrs_nnn),",
		"            where wavelength 'nnn' \u{2208} {412,443,469,488,531,547,555,645,667,678};",
		"  chl-a  -- chlorophyll-A concentration (CHL_chlor_a);",
		"  AOT      -- aerosol optical thickness, \u{03C4}, in one NIR or red band (RRS_aot_869);",
		"  \u{00C5}      -- angstrom coefficient (RRS_angstrom);",
		"  K490   -- diffuse attenuation coefficient at 490 nm (KD490_Kd_490);",
		"  PIC    -- calcite concentration or particulate inorganic carbon (PIC_pic);",
		"  POC    -- particulate organic carbon (POC_poc);",
		"  PAR    -- photosynthetically available radiation (PAR_par);",
		"  FLH    -- fluorescence line height (FLH_ipar, FLH_nflh);",
		"  IOPs   -- inherent optical properties (IOPs, numerous), which include",
		"            absorption and backscattering coefficients in the visible bands.",
		"  SST    -- sea surface temperature (NSST_sst, SST_sst, SST4_sst4);",
		"----------",
		"The MODIS SST products include 4-micron (nighttime only) and 11-micron (daytime and nighttime) SST.",
		"For the Level-3 products, each binned product contains multiple geophysical parameters,",
		"while the standard mapped image (SMI) products contain one parameter per granule.",
		"----------",
		"***** PBSsatellite cannot open the binned datasets, use the SMI products *****"
		)
		if (verbose)
			.flush.cat(paste0(mess,collapse="\n"),"\n\n")
		if(missing(years) || (missing(days) && missing(months)) || missing(res))
			stop("Supply 'years', 'days'/'months', and 'res' (resolution: '4km' or '9km') for this URL.")
		pname.choice = c("CHL_chlor_a", "RRS_aot_869", "RRS_angstrom", "KD490_Kd_490)", "PIC_pic", "POC_poc", "PAR_par", "FLH_ipar", "FLH_nflh",
			"NSST_sst", "SST_sst", "SST4_sst4",
			paste0("RRS_Rrs_",c(412,443,469,488,531,547,555,645,667,678)))
		if(missing(pname) || !is.element(pname, pname.choice))
			stop(paste0("Choose 'pname' from:\n", paste0(pname.choice,collapse="\n")))
		period.choice = c("DAY", "8D", "MO", "R32", "YR")
		if(missing(period) || !is.element(period, period.choice))
			stop(paste0("Choose 'period' from: '", paste0(period.choice,collapse="', '"), "'"))
		if (!is.element(period,c("DAY","MO")))
			stop("R code needs revision to get file name correct") ## Still needs work
		sensor.choice = c("Aquarius", "CZCS", "MERIS", "MODIS-Aqua", "MODIS-Terra", "OCTS", "SeaWiFS", "VIIRS-JPSS1", "VIIRS-SNPP")
		names(sensor.choice) = c("Q", "C", "M", "A", "T", "O", "S", "V1", "V2")
		if(missing(sensor) || !is.element(sensor, sensor.choice))
			stop(paste0("Choose 'sensor' from:\n", paste0(sensor.choice,collapse="\n")))
		cnt  = 0
		for(i in years) {
			if (period=="MO"){
				mos  = formatC(months, width=2, format="d", flag="0")
				j0   = as.Date(paste0(i,"-01-01"))
				days = julian(as.Date(paste0(i,"-",mos,"-01")), origin=j0) + 1
				emon = max(months)
				eday = julian(as.Date(paste0(ifelse(emon==12,i+1,i),"-",ifelse(emon==12,1,emon+1),"-01")), origin=j0) 
				nday = diff(c(days,eday+1)); names(nday) = days
#browser();return()
			}
			for(j in days) {
				jj    = formatC(j, width=3, format="d", flag="0")  ## pad with zeroes
				jjj   = as.character(j)
				if (period=="DAY"){
					fname = paste0("A", i, jj, ".L3m_", period, "_", pname, "_", res, ".nc")
				} else if (period=="MO"){
					jje = formatC(j + nday[jjj] - 1, width=3, format="d", flag="0")
					fname = paste0("A", i, jj, i, jje, ".L3m_", period, "_", pname, "_", res, ".nc")
				} else {
					stop("Period chosen is not supported")
				}
				fpath = paste0(save.dir, "/", fname)
				## oceandata has a specific download address for files without using subfolders:
				web   = paste(url, fname, sep="/")  ## cgi/getfile
				dweb  = paste(sub("/cgi/getfile","",url), sensor, "L3SMI", i, jj, sep="/")  ## display folder for listing

				if(!file.exists(fpath) || overwrite) {
					## If user does not have a curl installation
					if (missing(curlbin)) {
						try.out = try(download.file(web, destfile=fpath, method="auto", mode="wb"), silent=TRUE)
						if (class(try.out)!="try-error")
							cnt = cnt + 1
					}
					## If user has specified a curl installation bin path
					else {
						cmd   = paste0(curlbin,"/curl --list-only ", dweb)
						wdump = system (cmd, intern=T)
						if (any(grepl(fname,wdump))) {
#browser();return()
							download.file(web, destfile=fpath, method="auto", mode="wb")
							cnt = cnt + 1
						} ## end if exists
					} ## end if curlbin
				} ## end if write|overwrite
			} ## end j loop (days)
		} ## end i loop (years)
		n.possible <- length(years) * length(days)
	} ## end if oceandata

	## -----------------
	## NCEI SST data
	## -----------------
	if (grepl("ncei",url))
	{
		## Suggestion by Chris Rooper (DFO, PBS)
		## DOWNLOAD RELEVANT MONTHS AND YEARS OF DATA
		year<-seq(1982,2017,1)

		for(i in years) {
			for(j in months) {
				jj <- formatC(j, width=2, format="d", flag="0")
				urlname <- paste(url,"ersst.v5.", i, jj, ".nc",sep="")  ## no such files
				download.file(urlname, destfile=paste(save.dir,"/ersst.v5.", i, jj, ".nc", sep=""), mode="wb", quiet=FALSE)
			}
		}
	}

	## Check if all data files were downloaded
	n.exist    <- length(list.files(save.dir))
	cat(cnt, "files downloaded of", n.possible, "requested", "\n")
	cat(n.exist, "files exist in", save.dir, "\n")
	return("something")
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~download.ncdfData


## extractSlices------------------------2015-04-10
## Extract slices from an ncdfData object.
## ---------------------------------------------NL
extractSlices <- function(ncdfData, slices=NULL, dates=NULL, tlim=NULL)
{
	## validate ncdfData is a ncdfData object
	if (class(ncdfData) != "ncdfData") {
		stop("ncdfData must be of class ncdfData")
	}
	## extractSlices requires either slices, dates, or tlim argument
	## and only one to be specified
	## TRUE == 1, so if a value is != 1 user has provided incorrect arguments
	if (sum(c(!is.null(slices), !is.null(dates), !is.null(tlim))) != 1) {
		stop("You have not provided the correct number of arguments.")
	}
	if(!is.null(slices)){
		## make backup of slice names and attributes
		sliceNames <- names(ncdfData)

		## convert names to NULL, will replace later as slices are likely to change
		names(ncdfData) <- NULL
		attributesNcdfData <- attributes(ncdfData)

		## ensures numeric vector
		if(class(slices) != "numeric"){
			stop("slices must be a vector of type numeric")
		}
		## subset ncdfData to the proper slices
		ncdfData <- ncdfData[slices]

		## restore attributes and new names to ncdfdata
		attributes(ncdfData) <- attributesNcdfData
		names(ncdfData) <- sliceNames[slices]
	} else if(!is.null(dates)) {
		## make backup of slice names and attributes
		sliceNames <- names(ncdfData)

		## convert names to NULL, will replace later as slices are likely to change
		names(ncdfData) <- NULL
		attributesNcdfData <- attributes(ncdfData)

		## ensure user has provide date strings
		if(class(dates) != "character") {
			stop("dates must be a vector of type character")
		}
		## find the names of slices that are equal to the ones provided in dates
		## sort them to keep slices in order
		slicesOfDates <- sort(match(dates, sliceNames))
		ncdfData <- ncdfData[slicesOfDates]

		## restore attributes and correct slice names
		attributes(ncdfData) <- attributesNcdfData
		names(ncdfData) <- sliceNames[slicesOfDates]
	} else if(!is.null(tlim)) {
		ncdfData <- .ncdfDataTlimClip(ncdfData, tlim)
		## user has not provided either a slices, dates, or tlim argument 
	} else {
		stop("Must provide either slices or dates or tlim")
	}
	return(ncdfData)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~extractSlices


## extractTimeSeries--------------------2015-04-10
## Extract time series from an ncdfData object.
## ---------------------------------------------NL
extractTimeSeries <- function(ncdfData, xlim=NULL, ylim=NULL,
   polygons=NULL,functions=c("sum", "mean", "sd"),
   na.rm=TRUE, tlim=NULL, combine=1, by=NULL, 
   include.lowest=TRUE)
{
	## validate ncdfData is a ncdfData object
	if (class(ncdfData) != "ncdfData") {
		stop("ncdfData must be of class ncdfData")
	}
	## if combine is not a factor of length(ncdfData), remove excess slices, warn user
	if(combine > 1 && length(ncdfData)%%combine != 0){
		remainderOfSlices <- length(ncdfData)%%combine
		maxSlices <- length(ncdfData) - remainderOfSlices
		ncdfData <- extractSlices(ncdfData, slices=as.numeric(c(1:maxSlices)))
		warning("Combine not a factor of length(ncdfData), removing slices at end of ncdfData")
	}
	## if user has provided a tlim to clip ncdfData based on time
	if(!is.null(tlim)) {
		ncdfData <- .ncdfDataTlimClip(ncdfData, tlim)
	}
	
	if(!is.null(by)) {
		## backup ncdfData attributes
		ncdfDataAttributes <- attributes(ncdfData)

		## save sequence of indicated slices
		sliceSequence <- seq(from=1, to=length(ncdfData), by=by)

		## convert names to be consistent with new slices
		ncdfDataAttributes$names <- names(ncdfData)[sliceSequence]

		## convert ncdfData to have on the indicated slices 
		ncdfData <- ncdfData[sliceSequence]
		attributes(ncdfData) <- ncdfDataAttributes
	}
	## if user has specified a xlim argument
	if(!is.null(xlim)) {
		## function returns ncdfData with clipped xlim
		ncdfData <- .ncdfDataXClip(ncdfData, xlim, include.lowest)
	}
	## if user has specified a ylim argument 
	if(!is.null(ylim)) {
		## function returns ncdfData with clipped ylim
		ncdfData <- .ncdfDataYClip(ncdfData, ylim, include.lowest)
	}
	## if no polygon, convert xlim and ylim into a square polygon
	## is this an acceptable way to deal with this?
	## i was thinking it was the best way to take one singular/smaller path 
	## which we discussed as being beneficial in our meetings 
	if(is.null(polygons)) {
		xvect <- attr(ncdfData, "x")
		yvect <- attr(ncdfData, "y")
		x1 <- max(xvect)
		x2 <- min(xvect)
		y1 <- max(yvect)
		y2 <- min(yvect)

		## creates 1 polygon spanning all of xlim and ylim
		polygons <- data.frame(PID=c(rep(1, 4)), POS=c(1:4), X=c(x1, x1, x2, x2), Y=c(y2, y1, y1, y2))
	}
	## convert first slice to EventData
	## note: only need to clip one slice and then apply clipping result
	## to all existing slices
	ed <- to.EventData(ncdfData, 1)
	pointsInPoly <- findPolys(ed, polygons)
	if(is.null(pointsInPoly)) {
		stop("Your xlim and/or ylim has clipped out all of your polygons")
	}
	## numberSequenceInfo is a list that holds a numberSequence used to form variable 
	## length polygons into the same vector for each polygon when using combine. Also 
	## contains a vector of valid polygons that are being used, it is possible for polygons 
	## to be missing due to xlim and ylim arguments. If combine is 1, numberSequence will be 
	## an empty vector but the list will contain validPolyons.
	numberSequenceInfo <- .getNumberSequenceInfo(pointsInPoly, combine)

	## contains list of each polygon with its EID's/indexes
	polygonList <- split(pointsInPoly, pointsInPoly[, "PID"])

	dataInEachPoly <- lapply(ncdfData,".getPolysVectorList", polygonList)
	## dataInEachPoly is a list of slice dates with a list of each polygons data vector
	## dataInEachPoly structure
	## 1989-12-31: List of 2 (2 polygons)
	## $ 1: num [1:565] -1.8 -1.8 -1.8 -1.8 data vector of points in first PID
	## $ 2: num [1:400] 15.3 15.5 15.8 15.6 data vector of points in second PID
	## 1990-01-07: List of 2
	## ...

	if(combine != 1 && combine > 0) {
		## valuesVector contains all of the data for each lat/lon position in each polygon
		## for every slice in a single consecutive vector
		valuesVector <- unlist(dataInEachPoly)

		## do not need the names at this point, might save and assign later
		names(valuesVector) <- NULL

		## get the length of each slice with its data in each polygon
		lengthSliceValues <- length(valuesVector) / length(ncdfData)

		## get the length to split each combined slice, based on combine value
		lengthToSplit <- lengthSliceValues * combine

		## number of list values (slice) after the combine function
		numberOfLists <- length(valuesVector) / lengthToSplit

		## only works if combine goes in evenly, this might be a problem....
		## is there a way around this?
		## combine needs to be a factor of length(ncdfData), could always use extract slices
		## combinedSlicesList contains grouped slices (combine amount) with polygon value   

		combinedSlicesList <- split(valuesVector, rep(1:numberOfLists, each=lengthToSplit))

		numberSequence <- numberSequenceInfo$numberSequence

		## returns a list of data frames of combined slices with function results
		dataFramesList <- lapply(combinedSlicesList, function(combinedSlice, polyCount, combine, numberSequence)
		{
			## length of a singular combined slice
			lengthOfSlice <- length(combinedSlice) / combine

			## puts combinedSlice into a matrix
			dim(combinedSlice) <- c(lengthOfSlice, combine)
			## combinedSliceVector contains vector with ordered polygon data
			## e.g., if combine = 2 it would contain elements from
			## polygons 1,1,2,2
			combinedSliceVector <- unlist(t(combinedSlice))

			## splitting by numberSequence combines vectors of all of the same polygon and puts them
			## into a an order list starting with the first polygon and ending with the last in 
			## polygons
			combinedListOfPolygons <- split(combinedSliceVector, numberSequence)
			return(combinedListOfPolygons)
		}, length(unique(polygons$PID)), combine, numberSequence) 

		timeSeries <- lapply(dataFramesList,".getFunctionResults",functions, numberSequenceInfo$validPolys, na.rm)

		## gives time series inclusive date ranges
		names(timeSeries) <- paste(names(ncdfData)[seq(1, length(ncdfData), combine)],
			names(ncdfData)[c(seq(1, length(ncdfData), combine)[-1]-1, length(ncdfData))])
	} else {
		timeSeries <- lapply(dataInEachPoly,".getFunctionResults", functions, numberSequenceInfo$validPolys, na.rm)
		names(timeSeries) <- names(ncdfData)
	}
	return(timeSeries)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~extractTimeSeries


## listToDF-----------------------------2015-04-10
## Convert a list to a data.frame
## ---------------------------------------------NL
listToDF <- function(lst, newColumnName = "names")
{
	## validate our assumptions for the list of data frames
	## all data frames must have the same number of columns
	nc <- unique(unlist(lapply(lst, ncol)))
	if (length(nc) > 1) {
		stop("all data frames must have the same number of columns")
	}
	## all data frames must have the same column names
	colNameMatrix <- matrix(unlist(lapply(lst, "names")), ncol=nc, byrow=T)
	if (sum(!duplicated(colNameMatrix)) != 1) {
		stop("all data frames must have the same column names in the same order")
	}
	df <- do.call("rbind", lst)
	rownames(df) <- NULL

	## ensure that the new date column's name doesn't already exist
	if (is.element(newColumnName, names(df))) {
		stop("cannot add slice column; already exists; use sliceColumn to override")
	}
	## account for the possibility that a list element has a NULL value by
	## created a times vector the same length as names (init. to 0); whereever
	## nrows will only evaluate to an integer where it isn't NULL
	times <- numeric(length(names(lst)))
	times[!unlist(lapply(lst, "is.null"))] <- unlist(lapply(lst, "nrow"))
	df[[newColumnName]] <- rep(names(lst), times=times)

	colNames <- names(df)
	df <- df[, c(colNames[length(colNames)], colNames[-length(colNames)])]
	return (df)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~listToDF


## merge.ncdfData-----------------------2018-11-05
## Read a suite of netCDF binary files|objects,
## merge them, and create one ncdfData object in R.
## ---------------------------------------------RH
merge.ncdfData <- function (filenames, 
   dataVariable=1, convertMissingValues=FALSE,
   dataType=NULL, dataUnits=NULL, dataVname=NULL,
   xlim=NULL, ylim=NULL, tlim=NULL,
   x=NULL, y=NULL, time=NULL,
   Ux=NULL, Uy=NULL, Utimes=NULL, verbose=FALSE,
   read.ncdf=TRUE, write.ncdf=FALSE, outname)
{
	on.exit(gc(verbose=FALSE))
	flist = list()
	nfile = length(filenames)
	if (read.ncdf) {
		for (f in 1:nfile) {
			ff = filenames[f]
#browser();return()
			if (!file.exists(ff)) next
			if (!is.null(Utimes)) Utime = (rep(Utimes,nfile))[f]
			else                  Utime = NULL
			flist[[ff]] = read.ncdfData (ff, 
				dataVariable = dataVariable,
				convertMissingValues = convertMissingValues,
				dataType = dataType,
				dataUnits = dataUnits,
				xlim=xlim, ylim=ylim, tlim=tlim,
				x=x, y=y, time=time,
				Ux=Ux, Uy=Uy, Utime=Utime, verbose=verbose
			)
		}
	} else {  ## just populate flist with existing ncdfData objects
		for (f in 1:nfile) {
			ff = filenames[f]
			if (!exists(ff) || class(get(ff))!="ncdfData") next
			flist[[ff]] = get(ff)
		}
	}
	if (length(flist)==0) stop("No valid ncdf data files|objects were collected")
	
	## Assume the first ncdfData object in the series represents subsequent ones
	f1.atts   = attributes(flist[[1]])
#browser();return()
	f1.atts   = f1.atts[setdiff(names(f1.atts),"names")]

	## If specified by the user, overide existing descriptors
	if (!is.null(dataType)) {
		f1.atts$dataType <- dataType
	}
	if (!is.null(dataUnits)) {
		f1.atts$dataUnits <- dataUnits
	}
	if (!is.null(dataVname)) {
		f1.atts$dataVname = dataVname
	}

	## Build a merged ncdfData set
	ncdata = list()
	for (f in 1:length(flist)) {
		ff    =  names(flist[[f]])
		fdata = flist[[f]][[ff]]
		ncdata[[ff]] = fdata
	}
	## Add attributes to ncdfData object
	nc.atts = attributes(ncdata)
	attributes(ncdata) = c(nc.atts, f1.atts)

	if (write.ncdf) {
		## Create a netCDF file from merged files
		if (missing(outname))
			outname = "ncdf.merged.nc"
		Tvals=names(ncdata)
		tvals = as.numeric(as.Date(Tvals) - as.Date(Tvals[1]))
		tmess = paste0("days since ",Tvals[1])

		## Construct zfld array
		xvect = f1.atts$x; yvect = f1.atts$y
		nx = length(xvect); ny = length(yvect); nt = length(tvals)
		zmat  = array(NA, dim=c(nx,ny,nt), dimnames=list(x=xvect,y=yvect,t=Tvals))
		for (i in Tvals) {
			imat = ncdata[[i]][["data"]]
			dimnames(imat) = list(x=xvect,y=yvect)
			zmat[dimnames(imat)[[1]], dimnames(imat)[[2]], i] = imat
		}
		zlist = list()
		zlist[[f1.atts$dataVname]] = zmat
#browser();return()
		create.ncdfData(filename=outname, xvals=xvect, yvals=yvect, tvals=tvals, tmess=tmess, dataType=f1.atts$dataType, dataUnits=f1.atts$dataUnits, zlist=zlist)
	}
	return(ncdata)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~merge.ncdfData


## plot.ncdfData------------------------2018-11-08
## Plot an ncdfData slice
## ------------------------------------------NL|RH
plot.ncdfData <- function(x, slice, layer="data",
   xlim=NULL, ylim=NULL, style=c("image","contour"), projection="LL",
   tck=-0.014, tckMinor=0.5*tck, showXY=FALSE, ...)
{
	if (missing(slice)) {
		warning("slice argument missing; assuming first slice")
		slice <- 1
	}
	style <- match.arg(style)
	slice <- .getSlice(x, slice)

	xvect <- attributes(x)$x
	yvect <- attributes(x)$y

	z <- slice[[layer]]
	if (is.null(z)) {
		stop(paste0("The layer '", layer, "' does not exist in the ncdfData object."))
	}

	## image(...) and contour(...) both require both x and y to be in
	## *increasing* order;  our standardized representation has y in
	## *decreasing* order ==> fix it
	yvect <- rev(yvect)
	z <- z[, ncol(z):1, drop=FALSE]

	## plot the whole range by default
	if (is.null(xlim)) {
		xlim <- range(xvect)
		if (diff(xlim)==0)
			xlim = xlim + ( ifelse(is.null(attributes(x)$xBy),1,attributes(x)$xBy/2) * c(-1,1) )
	}
	if (is.null(ylim)) {
		ylim <- range(yvect)
		if (diff(ylim)==0)
			ylim = ylim + ( ifelse(is.null(attributes(x)$yBy),1,attributes(x)$yBy/2) * c(-1,1) )
	}

	## some argument names conflict between plotMap and image/contour;
	## remove those arguments from '...' and pass them only into
	## image/contour
	args <- list(...)
	args <- args[!is.element(names(args), c("col", "lty", "lwd","zlim"))]

	## add required arguments for plotMap
	args <- c(list(polys=NULL), args) # otherwise, removes polys from list
	args$xlim <- xlim
	args$ylim <- ylim
	args$projection <- projection
	args$tck <- tck
	args$tckMinor <- tckMinor
	do.call("plotMap", args)

	## pass all arguments from '..' to contour/image
	args <- list(...)
	## add required arguments for contour/image
	args$x <- if (length(xvect)>1) xvect else xvect + ( ifelse(is.null(attributes(x)$xBy),1,attributes(x)$xBy/2) * c(-1,1) )
	args$y <- if (length(yvect)>1) yvect else yvect + ( ifelse(is.null(attributes(x)$yBy),1,attributes(x)$yBy/2) * c(-1,1) )
	args$z <- z
	args$add <- TRUE

	do.call(style, args)

	if (showXY) { ## (RH 181022)
		## Note: ncdfData made from irregularly spaced data will have differing grid cell sizes.
		xymat = args$z; rownames(xymat) = args$x; colnames(xymat) = args$y
		xyuse = which(!is.na(xymat), arr.ind=TRUE)
		XY    = apply(cbind(X=rownames(xyuse), Y=colnames(xymat)[xyuse[,"col"]]), 2, as.numeric)
		if (is.null(dim(XY)))
			XY = matrix(XY,nrow=1); colnames(XY)=c("X","Y")
		xypts = as.EventData(data.frame(EID=1:nrow(XY), XY), projection="LL")
		addPoints(xypts,col="green",pch=20)
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plot.ncdfData


## print.ncdfData-----------------------2018-10-31
## Print a summary of the ncdfData object.
## ------------------------------------------NL|RH
print.ncdfData <- function(x, ...)
{
	## extract some data to simplify printing below
	xCoords <- attributes(x)$x
	yCoords <- attributes(x)$y

	xRng <- sort(range(xCoords))
	yRng <- sort(range(yCoords))

	## added called to round(...) because the difference may not be
	## as even as we might expect; in chlonasaclim.nc from
	## http://research.jisao.washington.edu/data_sets/seawifs/
	## observed the following unique differences:
	## 1.200000 1.200000 1.200000 1.200000 1.200001 1.199999 1.200003 1.199997
	## 1.200005 1.200012 1.199982
	if (length(xCoords)==1)
		xBy = ifelse(is.null(attributes(x)$xBy),0,attributes(x)$xBy)
	else
		xBy <- abs(unique(round(diff(xCoords), 3)))
	if (length(yCoords)==1)
		yBy = ifelse(is.null(attributes(x)$yBy),0,attributes(x)$yBy)
	else
		yBy <- abs(unique(round(diff(yCoords), 3)))

	cat(sprintf("NCDF data\n"))
	cat(sprintf("\tData type:     %s\n", attributes(x)$dataType))
	cat(sprintf("\tData units:    %s\n", attributes(x)$dataUnits))
	cat(sprintf("\tVariable name: %s\n", attributes(x)$dataVname))

	cat(sprintf("\nSlices:\n"))
	cat(sprintf("\tCount: %d\n", length(attributes(x)$names)))
	cat(sprintf("\tFirst: %s\n", head(attributes(x)$names, 1)))
	cat(sprintf("\tLast:  %s\n", tail(attributes(x)$names, 1)))
	
	cat(sprintf("\nSlice data:\n"))
	cat(sprintf("\tX: %8.3f to %8.3f by %6.3f\n", xRng[1], xRng[2], xBy))
	cat(sprintf("\tY: %8.3f to %8.3f by %6.3f\n", yRng[1], yRng[2], yBy))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~print.ncdfData


## read.ncdfData------------------------2018-11-05
## Read a netCDF binary file and create an
## ncdfData object in R
## ------------------------------------------NL|RH
read.ncdfData <- function (filename, 
   dataVariable=1, convertMissingValues=FALSE, 
   dataType=NULL, dataUnits=NULL, dataVname=NULL,
   xlim=NULL, ylim=NULL, tlim=NULL,
   x=NULL, y=NULL, time=NULL,
   Ux=NULL, Uy=NULL, Utime=NULL, verbose=FALSE)
{
	if (!file.exists(filename))
		stop(paste0(filename, " does not exist"))
	nc <- nc_open(filename, readunlim=FALSE)
	on.exit(nc_close(nc))

	.flush.cat = function (...) {
		cat(...); flush.console(); invisible() }

	## Let user know what variables are available
	if (verbose) {
		mess = paste0("Verbosity:\n---------\nFile '", filename, "'\ncontains ", nc$nvars, " variables: '", paste0(names(nc$var), collapse="', '"), "'")
		.flush.cat(mess, "\n\n")
	}
	## Get the data object specified by the user
	data    = nc$var[[dataVariable]]
	dnams   = sapply(nc$var[[dataVariable]]$dim,function(x){x$name})
	ddims   = sapply(nc$var[[dataVariable]]$dim,function(x){length(x$vals)})
	if (verbose) {
		mess = paste0("'",dnams,"'[",ddims,"]",collapse=", ")
		mess = paste0("Variable '", names(nc$var)[dataVariable], "' contains ", length(dnams), " dimensions: ", mess)
		.flush.cat(mess, "\n\n")
	}
	## If unspecified by the user (NULL), set sensible defaults
	if (is.null(dataType)) {
		dataType <- data$longname
	}
	if (is.null(dataUnits)) {
		dataUnits <- data$units
	}
	if (is.null(dataVname)) {
		dataVname = names(nc$var)[dataVariable]
	}
	if (is.null(x)) {
		z = grep("^X|^[Ll]on",dnams)
		if (length(z)==0) x = "lon"
		else x = nc$var[[dataVariable]]$dim[[z]]$name
	}
	if (is.null(y)) {
		z = grep("^Y|^[Ll]at",dnams)
		if (length(z)==0) y = "lat"
		else y = nc$var[[dataVariable]]$dim[[z]]$name
	}
	if (is.null(time)) {
		z = grep("^[Tt]ime|^[Dd]ate",dnams)
		if (length(z)==0) time = "time"
		else time = nc$var[[dataVariable]]$dim[[z]]$name
	}
	## Anytime the user did not specify a dimension's unit,
	## attempt to extract it automatically
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
	## If units not found, we also could not find the dimension
	if (is.null(Ux)) stop("x dimension not found; specify using 'x'")
	if (is.null(Uy)) stop("y dimension not found; specify using 'y'")
	if (is.null(Utime)) stop("time dimension not found; specify using 'time'")

	## Store x and y vectors containing the spatial representation of data points
	xvect <- as.vector(ncvar_get(nc,x))
	yvect <- as.vector(ncvar_get(nc,y))

	## Store length of x and y vectors
	xlen <- length(xvect)
	ylen <- length(yvect)

	## Split up time units into a vector
	## e.g., time units might be "days since 1900-01-01 00:00:00"
	timeV <- strsplit(Utime, " ")[[1]]

	if (time %in% dnams) {
		## Get the actual time stamps from nc file
		raw.names <- ncvar_get(nc, time)
		## Convert time stamps to strings in the format
		## "YYYY-MM-DD HH:MM:SS"
		names <- .getSliceNames(raw.names, unit=timeV[1], origin=paste(timeV[3], timeV[4]))
	} else {
		names = raw.names = rev(timeV)[1]
	}
	if (verbose) {
		mess = paste0("Utime: ", Utime, "\nraw.names: ", paste(raw.names,collapse=" "), "\nnames:\n", paste(names,collapse="\n"), collapse="")
		.flush.cat(mess, "\n\n")
	}
	## User has specified a tlim argument, need to limit the amount of slices
	if (!is.null(tlim)) {
		if(!is.character(tlim)){
			stop("'tlim' must be of type character please specify a date string e.g., '2000-01-01'")
		}
		if (length(tlim) != 2) {
			stop("'tlim' must be a vector of two strings, each formatted as \"YYYY-MM-DD HH:MM:SS\"")
		}
		## Convert date strings into POSIXlt objects for date comparison
		names <- as.POSIXlt(names, tz="GMT", origin=paste(timeV[3], timeV[4]))

		## Keep only the dates that fall in tlim
		datesInTLIM <- names >= as.POSIXlt(tlim[1], tz="GMT") &
			names <= as.POSIXlt(tlim[2], tz="GMT")

		## Get locations of the data sheets in the data object
		rangeT <- which(datesInTLIM)

		## Use locations of first and last value in rangeT to create of offset into the data object
		## this is done to prevent unnecessary processing of data that is not in tlim
		dataset <- ncvar_get(nc, data,
							 start=c(1, 1, min(rangeT)),
							 count=c(xlen, ylen, max(rangeT)-min(rangeT)+1))

		## Convert POSIXlt objects back to strings for name assignment
		names <- as.character(names[datesInTLIM])

		## No tlim specified, read in the full contents of data object
	} else {
		dataset <- ncvar_get(nc, data)
	}
	## Ensure that dataset has 2 or 3 dimensions -- assumed by the code below
	if (length(dim(dataset)) < 2 || length(dim(dataset)) > 3) {
		stop("Unsupported data variable: contains either less than 2 or more than 3 dimensions")
	}
	## Add a third dimension to two-dimensional objects to simplify processing below
	if (length(dim(dataset)) == 2) {
		dim(dataset) <- c(dim(dataset), 1)
	}
	
	## Standardize the order for x and y and reorder columns/rows as necessary:
	## - xvect should be in *increasing* order
	## - yvect should be in *decreasing* order
	## when indexing dataset, dataset[1, 1] is thus the top-left of the map
	## meaning min(x), max(y); this representation is consistent with the majority
	## of encountered NetCDF data sets
	if (head(xvect, 1) > tail(xvect, 1)) {
		xvect <- rev(xvect)
		dataset <- dataset[dim(dataset)[1]:1, , , drop=FALSE]
	}
	if (head(yvect, 1) < tail(yvect, 1)) {
		yvect <- rev(yvect)
		dataset <- dataset[, dim(dataset)[2]:1, , drop=FALSE]
	}
	## Clip the data to keep only data within xlim/ylim
	if (!is.null(xlim)) {
		xsub <- xvect >= xlim[1] & xvect <= xlim[2]  ## Just use logical vector (RH)
		xlen <- sum(xsub)
		if (xlen<1) stop ("Subsetting using `xlim' yields zero records")
		xvect = xvect[xsub]
		dataset <- dataset[xsub, , , drop=FALSE]
	}
	if (!is.null(ylim)) {
		ysub <- yvect >= ylim[1] & yvect <= ylim[2]  ## Just use logical vector (RH)
		ylen <- sum(ysub)
		if (ylen<1) stop ("Subsetting using `ylim' yields zero records")
		yvect = yvect[ysub]
		dataset <- dataset[ ,ysub, , drop=FALSE]
	}
#browser();return()
	## Convert NAs back to native missval if specified
	## nc_open automatically converts missval to NA
	## this gives the users the choice to switch back if they need
	if (convertMissingValues == TRUE) {
		missval <- data$missval
		missingLocations <- which(is.na(dataset))
		dataset[missingLocations] <- data$missval
	} else {
		missval <- NA
	}
	## Calculate size of each "data sheet" in the nc file
	size <- (xlen * ylen)

	## Calculate number of data sheets in the dataset
	length <- length(dataset) / size

	## Convert dataset into a list of vectors
	ncdfData <- split(dataset, rep(1:length, each=size))
#browser(); return()

	## Put dimension on each list vector to convert to a list of matrices;
	## names each matrix as data
	ncdfData <- lapply(ncdfData, ".createDim", xlen, ylen)

	## Apply name changes to the list
	names(ncdfData) <- names

	## add attributes to ncdfData object
	attr(ncdfData,"dataType")  <- dataType
	attr(ncdfData,"dataUnits") <- dataUnits
	attr(ncdfData,"dataVname") <- dataVname
	attr(ncdfData,"x")         <- xvect
	attr(ncdfData,"y")         <- yvect
	attr(ncdfData,"missval")   <- missval

	#nc_close(nc) ## uses `on.exit' function to close the nc file

	## give class name to ncdfData object
	class(ncdfData) <- "ncdfData"
	return(ncdfData)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~read.ncdfData


## removeAnomalousValues----------------2015-04-10
## Remove anomalous values from ncdfData
## ---------------------------------------------NL
removeAnomalousValues <- function(ncdfData, zlim)
{
	if (length(zlim) != 2) {
		stop("zlim must be a range of two values; use NA to omit a minimum/maximum")
	}
	## make backup of ncdfData attributes
	ncdfDataAttributes <- attributes(ncdfData)

	ncdfData <- lapply(ncdfData, function(slice, zlim) {
		if (!is.na(zlim[1])){
			slice$data[slice$data < zlim[1]] <- NA
		}
		if (!is.na(zlim[2])){
			slice$data[slice$data > zlim[2]] <- NA
		}
		return(slice)
	}, zlim)

	## restore attributes
	attributes(ncdfData) <- ncdfDataAttributes
	return(ncdfData)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~removeAnomalousValues


## scaleRegion--------------------------2015-04-10
##  Scale ncdfData slices to a new resolution 
##  based on a scale factor.
## ---------------------------------------------NL
scaleRegion <- function(ncdfData, scaleFactor, fun="drop",
   placement="topleft", includeErrorMatrix=FALSE,
   includeMissMatrix=FALSE, remainder="crop", na.rm=TRUE)
{
	## validate ncdfData is a ncdfData object
	if (class(ncdfData) != "ncdfData") {
		stop("ncdfData must be of class ncdfData")
	}
	## backup ncdfData attributes
	attributesNcdfData <- attributes(ncdfData)

	## user provided a positive scale factor. The ncdfData object slices will be
	## scaled up based on the scale factor. For PBSsatellite version one repeat
	## method will be used.  The number of points in each slice will be
	## increased by scaleFactor:1 ratio.
	if(scaleFactor > 0) {
		if(fun != "repeat"){
			warning("Version 1 always uses repeat method for scaling up.")
		}
		## create the scaled slices
		ncdfData <- lapply(ncdfData, ".createScaledUpSlices", scaleFactor)

		## restore the attributes
		attributes(ncdfData) <- attributesNcdfData

		## change x and y attributes to match the new ncdfData object
		x <- attr(ncdfData, "x")
		y <- attr(ncdfData, "y")

		## create vectors of new x and y variables
		x <- seq(from=min(x), to=max(x), length=length(x) * scaleFactor)

		## ensure y values are in decreasing order
		y <- sort(seq(from=min(y), to=max(y), length=length(y) * scaleFactor), decreasing=TRUE)

		## assign new x and y attributes to ncdfData object
		attr(ncdfData, "x") <- x
		attr(ncdfData, "y") <- y

		## user has provided a negative scale factor. The ncdfData
		## object slices will be scaled down  based on the method the user has
		## specified. The number of points each slice will be reduced by
		## 1/scaleFactor^2.
	} else {
		## retrieve both x and y ncdfData attributes
		x <- attr(ncdfData, "x")
		y <- attr(ncdfData, "y")

		## if scale factor is not a factor of length(x)
		## we need to add rows of NA to make it a factor
		if(length(x) %% scaleFactor != 0) {
			rowsToChange <- length(x) %% scaleFactor

			## rows to change can be either positive or negative
			## if positive add rows (fill) if negative removing rows (crop)
			if(remainder == "crop"){
				rowsToChange <- rowsToChange - abs(scaleFactor)
			}
			## add or remove rows to ncdfData slices
			ncdfData <- lapply(ncdfData, function(slice, rowsToChange)
			{
				if(rowsToChange > 0) {
					## add rowsToChange number of rows of NAs values to slices
					slice$data <- rbind(slice$data, matrix(NA, nrow=rowsToChange, ncol=ncol(slice$data)))
				} else {
					## remove rowsToChange number of rows to slices
					slice$data <- slice$data[1:(nrow(slice$data) - abs(rowsToChange)), 1:ncol(slice$data)]
				}
				return(slice)
			}, rowsToChange)
		}
		## if scale factor is not a factor of length(y)
		## we need to add cols of NA to make it a factor
		if(length(y) %% scaleFactor != 0) {
			colsToChange <- length(y) %% scaleFactor

			## cols to change can be either positive or negative
			## if positive add cols (fill) if negative removing cols (crop)
			if(remainder == "crop"){
				colsToChange <- colsToChange - abs(scaleFactor)
			}
			## add or remove cols to ncdfData slices
			ncdfData <- lapply(ncdfData, function(slice, colsToChange)
			{
				if(colsToChange > 0){
					slice$data <- cbind(slice$data, matrix(NA, nrow=nrow(slice$data), ncol=colsToChange))
				} else {
					slice$data <- slice$data[1:nrow(slice$data), 1:(ncol(slice$data) - abs(colsToChange))]
				}
				return(slice)
			}, colsToChange)
		}
		## restore ncdfData attributes
		attributes(ncdfData) <- attributesNcdfData

		## if else statement stops new attributes from including clipped regions
		## which would step outside of the world boundary
		if(remainder == "crop") {
			xlen <- nrow(ncdfData[[1]]$data)
			ylen <- ncol(ncdfData[[1]]$data)
		} else {
			xlen <- length(attr(ncdfData, "x"))
			ylen <- length(attr(ncdfData, "y"))
		}
		## update ncdfData x and y attributes to have new values due to scaling down
		attr(ncdfData, "x") <- x[seq(from=1, to=xlen, by=abs(scaleFactor))]
		attr(ncdfData, "y") <- y[seq(from=1, to=ylen, by=abs(scaleFactor))]

		## returns a scaleVector used to split ncdfData slices in scaleFactor^2
		## size chunks. These chunks are then sent to fun as provided the user
		## for a scaling method. The return value from the function then
		## replaces the points in the chunk.
		scaleVector <- .createScaleSplitVector(ncol(ncdfData[[1]]$data),
			nrow(ncdfData[[1]]$data), length(ncdfData[[1]]$data), abs(scaleFactor))

		## backup ncdfData attributes
		attributesNcdfData <- attributes(ncdfData)

		## assuming the user specified includeErrorMatrix and includeMissMatrix,
		## the resulting structure (scaledDownPoints) will look as follows:
		## $`1989-12-31`$`1`
		## $`1989-12-31`$`1`$data
		##[1] -1.8
		## $`1989-12-31`$`1`$error
		##[1] 0
		## $`1989-12-31`$`1`$miss
		##[1] 75
		## the `1` above is for the 4 (or 16 or ...) data points that were
		## collapsed into a single point
		scaledDownPoints <- lapply(ncdfData, ".createScaledDownSlices",
			scaleVector, scaleFactor, fun, na.rm, includeErrorMatrix, includeMissMatrix)

		## get the number of layers that will be in the new data slice, there
		## is always at least 1 (data) layer, and 2 optional error and/or
		## missing data layers
		numberOfLayers <- length(names(scaledDownPoints[[1]][[1]]))

		## layerSplitVector used to create either 1-3 layers in each ncdfData slice
		layerSplitVector <- rep(1:numberOfLayers, times=length(scaledDownPoints[[1]]))

		## apply regionSplitVector to each list element in scaledDownPoints
		## to create 1-3 layers and creates a new ncdfData object containing vectors.
		ncdfData <- lapply(scaledDownPoints, function(slice, layerSplitVector, layerNames, ncol, nrow)
		{
			## create a large vector containing (assuming all 3 layers)
			## point1:data,error,missing,point2:data,error,missing values
			slice <- unlist(slice)
			names(slice) <- NULL
			## split vectors into proper layer arrangement
			slice <- split(slice, layerSplitVector)
			## names the slices
			names(slice) <- layerNames

			## apply dimensions on slice to create matrix layout
			for(layerName in names(slice)){
				dim(slice[[layerName]]) <- c(ncol, nrow)
			}
			return(slice)
		}, layerSplitVector, names(scaledDownPoints[[1]][[1]]), length(attr(ncdfData,"x")), length(attr(ncdfData,"y")))

		## restore ncdfData attributes
		attributes(ncdfData) <- attributesNcdfData

		## change the attributes of ncdfData x and y to be shifted to the centre of each region
		if(placement == "centre") {
			if(remainder == "fill") {
				warning("x/y attrs may extend world boundary due to centre placement")
			}
			x <- attr(ncdfData, "x")
			y <- attr(ncdfData, "y")

			## find distance to centre of a point
			xDist <- abs(x[1] - x[2]) / 2
			yDist <- abs(y[1] - y[2]) / 2

			## change ncdfData attributes to be in the centre of scale down regions
			attr(ncdfData, "x") <- attr(ncdfData, "x") + xDist
			attr(ncdfData, "y") <- attr(ncdfData, "y") - yDist
		}
	}
	return(ncdfData)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~scaleRegion


## to.EventData-------------------------2015-04-10
## Convert an ncdfData object to an EventData object.
## ------------------------------------------NL|NB
to.EventData <- function(ncdfData, slice)
{
	## validate that ncdfData is a ncdfData object
	if(class(ncdfData) != "ncdfData"){
		stop("ncdfData must be of class ncdfData")
	}
	slice <- .getSlice(ncdfData, slice)
	
	## get lon and lat attributes from ncdfData
	lon <- attr(ncdfData,"x")
	lat <- attr(ncdfData,"y")

	## vector length for one week of data
	size <- length(lon) * length(lat)

	## build the data frame
	df <- data.frame(EID=1:size)
	df$X <- rep(lon, times=length(lat))
	df$Y <- rep(lat, each=length(lon))

	## if layerName left NULL use "DATA"

	## extract *all* of the layers -- making each layer into a column (column name = layer name)
	## for layerName in names(slices)...
	##   pull the layer's data out into a new column
	for (layerName in names(slice)) {
		df[[layerName]] <- as.vector(slice[[layerName]])
	}
	df <- as.EventData(df)
	return(df)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~to.EventData


## write.ncdfData-----------------------2018-11-09
## Write slices of ncdfData to CSV files.
## ---------------------------------------------RH
write.ncdfData = function(ncdfData, layer="data", slice, 
   file.prefix, save.dir="./ncdfData")
{
	if (class(ncdfData)!="ncdfData")
		stop("Supply a valid 'ncdfData' class object")
	if (missing(slice))
		slice = 1
	else if (is.character(slice) && all(slice=="all"))
		slice = 1:length(ncdfData)

	.getUdig = function(x) {
		dig = 0
		while (length(unique(round(x,dig)))!=length(x))
			dig = dig + 1
		return(dig)
	}
	xvals = attributes(ncdfData)$x
	xvec  = round(xvals, .getUdig(xvals) + 1)
	yvals = attributes(ncdfData)$y
	yvec  = round(yvals, .getUdig(yvals) + 1)
	snams = names(ncdfData)

	if (missing(file.prefix))
		file.prefix =  gsub(" ","_",attributes(ncdfData)$dataVname)

	## Get rid of trailing delimiters (if they are specified)
	save.dir = sub("/$","",save.dir)
	if(!dir.exists(save.dir))
		dir.create(save.dir, recursive = TRUE)

	for (i in slice) {
		if (is.numeric(i))
			ii = snams[i]
		else {
			if (!is.element(i, snams)) next
			ii = i
		}
		inam = paste0(save.dir,"/",file.prefix,"_",ii)
		imat = ncdfData[[i]][[layer]]
		rownames(imat) = xvec
		colnames(imat) = yvec
		write.csv(imat, paste0(inam,".csv"))
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~write.ncdfData


##-----Supplementary hidden functions-----

## .createDim---------------------------2015-04-10
## Applies a dimension to a list of vectors, puts the new matrix
## into a list and names the matrix data for the data layer
## ---------------------------------------------NL
.createDim <- function(vectors, xlen, ylen)
{
	## at this point, one would index the first dimension with
	## an x value and the second with a y value -- in the
	## order you would expect for indexing, i.e., (x, y)
	dim(vectors) <- c(xlen, ylen)
	## put matrix (vector) in a list, name matrix data
	vectors <- list(data=vectors)
	return(vectors)
}


## .createScaledDownSlices--------------2015-04-10
##  Helper function to scale down slices.
## ---------------------------------------------NL
.createScaledDownSlices <- function(slice, scaleVector, scaleFactor, fun,
   na.rm, includeErrorMatrix, includeMissMatrix)
{
	## creates a list of vectors containing each region of scaleFactor^2 points
	sliceRegionList <- split(slice$data, scaleVector)

	## apply functions and create additional layers if necessary
	scaledDownPointsList <- lapply(sliceRegionList, function(regionData, fun,
			vectorLength, na.rm, includeErrorMatrix, includeMissMatrix)
	{
		## user has provided a function other than the default drop
		if(fun != "drop"){
			## omitNA[1] == NA if all points in region data are NA
			omitNA <- na.omit(regionData)
			if(!is.na(omitNA[1])) {
				## remove NA values from region data
				if(na.rm) {
					cmd <- paste("newData <- ", fun[1], "(na.omit(regionData))", sep="")
					expression <- parse(text=cmd)
					eval(expression)
					## do not remove NA values from region data
				} else {
					cmd <- paste("newData <- ", fun[1], "(regionData)", sep="")
					expression <- parse(text=cmd)
					eval(expression)
				}
				## all values are NA
			} else {
				newData <- NA
			}
			## function is drop, keep the first point in regionData
		} else {
			newData <- regionData[1]
		}
		## regionList is a list of data points after function
		regionList <- list(data=newData)

		if(includeErrorMatrix) {
			## still looking into options for best error calculation
			## ask Rowan/Lyse how to calculate the error
			if(!is.na(newData)){
				error <- abs(regionData[1] - newData)
				percentError <- (error / abs(regionData[1])) * 100
				## if no value was calculated due to all NA's NA will be the percent error
			} else {
				percentError <- NA
			}
			regionList[["error"]] <- percentError
		}
		if(includeMissMatrix){
			## is.na includes nan
			## subtract nan count from na to get number of NA (missing data)
			nanCount <- sum(is.nan(regionData))
			naAndnanCount <- sum(is.na(regionData))
			## left with NA value count
			miss <- naAndnanCount - nanCount
			## calculate missing percentage
			miss <- (miss/vectorLength) * 100
			regionList[["miss"]] <- miss
		}
		return(regionList)
	}, fun, scaleFactor*scaleFactor, na.rm, includeErrorMatrix, includeMissMatrix)
	return(scaledDownPointsList)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~.createScaledDownSlices


.createScaledUpSlices <- function(slice, scaleFactor)
{
	## repeat a slice scaleFactor times and returns a vector of
	## scaleFactor combined slices
	scaledVector <- rep(slice$data, each=scaleFactor)
	## create a split factor used to split scaledVector into a
	## list of columns
	splitFactor <- rep(1:ncol(slice$data),
		   each=nrow(slice$data) * scaleFactor)
	## list of columns with slice data
	splitList <- split(scaledVector, splitFactor)
	## repeats the list scaleFactor times
	splitList <- rep(splitList, each=scaleFactor)
	## unlist returns a vector of a larger new slice with repeated values
	newSlice <- unlist(splitList)
	## remove undesired names
	names(newSlice) <- NULL
	## add new dimension to convert vector in to a matrix
	dim(newSlice) <- scaleFactor * dim(slice$data)
	return(list(data=newSlice))
}


## .createScaleSplitVector--------------2015-04-10
## Creates a vector of indices that will be used to extract data from
## a ncdfData slice. In order to scale down a ncdfData object data
## is needed from specified indices to properly scale down the object
## ---------------------------------------------NL
.createScaleSplitVector <- function(ncol, nrow, sizeSlice, scaleFactor)
{
	## get length of slices to be created
	newSliceSize <- sizeSlice / (scaleFactor * scaleFactor)

	## create matrix used for splitFactor
	m <- matrix(1:newSliceSize, ncol=ncol/scaleFactor)

	## repeat matrix scaleFactor to create a matrix 1/scaleFactor size of
	## original slice
	v <- rep(m, each=scaleFactor)

	## create a splitFactor used to create lists of index locations used for
	## scaling region
	splitFactor <- rep(1:ncol(m), each=nrow(m) * scaleFactor)

	## create a list of indexes used for acquiring correct slice data to be used
	## by the fun indicated to the user to accurately scale down a region
	vList <- split(v, splitFactor)

	## repeat the list of vectors by scaleRegion to accurately include all
	## indices in fun calculation, vList contains the same amount of data
	## indices as an original slice
	vList <- rep(vList, each=scaleFactor)

	## unlist vList to create a vector that will split a slice into the correct
	## locations used to be sent to a scaling fun to accurately scale down a region
	splitVector <- unlist(vList)
	names(splitVector) <- NULL

	return(splitVector)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~.createScaleSplitVector


## .findRC------------------------------2018-11-08
## Return number of rows and columns for plotting
## multi-panel figures given number of figures (nf)
## to fit on one page.
## Similar to function PBSmodelling::.findSquare
## ---------------------------------------------RH
.findRC = function (nf, orient="landscape") 
{
	sqn = sqrt(nf)
	m = ceiling(sqn)
	n = ceiling(nf/m)
	if (orient=="landscape")
		return(c(n, m))
	else
		return(c(m, n))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.findSquare


## .getClippedMatrix--------------------2015-04-10
## returns a mask of TRUE and FALSE values of the minimum matrix size
## to fit the full size of polygon.
## TRUE values are for values that are inside the polygon
## FALSE values are for values that are outside of the polygon.
## This result can be used to clip all slices in a ncdfData object
## for the given polygon.
## ---------------------------------------------NL
.getClippedMatrix <- function(eventData, xvect, yvect, polygons)
{
	## findPloys returns a LocationSet of polygons that are located inside
	## of polygon.
	pointsInPoly <- findPolys(eventData, polygons)
	
	if(is.null(pointsInPoly)){
		stop("You have clipped out all of your data")
	}
	## EID's are indexes of elements that belong inside the polygon
	indexesInPoly <- eventData[pointsInPoly[["EID"]], ]

	## creates a matrix the same size as slice with all FALSE values
	sliceMask <- matrix(FALSE, ncol=length(yvect), nrow=length(xvect))
	## replace locations inside of the polygon with TRUE
	sliceMask[indexesInPoly[["EID"]]] <- TRUE

	## find the dimensions of the smallest matrix that will fit polygon in its entirety
	newXVect <- xvect[xvect <= max(indexesInPoly["X"]) & xvect >=  min(indexesInPoly["X"])]
	newYVect <- yvect[yvect <= max(indexesInPoly["Y"]) & yvect >=  min(indexesInPoly["Y"])]

	## find the index locations of the minimum X value and maximum X value;
	## these values will be used to strip off unneeded columns of X values
	minX <- which(xvect==newXVect[1])
	maxX <- which(xvect==newXVect[length(newXVect)])

	## find the index locations of the minimum Y value and maximum Y value;
	## these values will be used to strip off unneeded rows of Y values
	minY <- which(yvect==newYVect[1])
	maxY <- which(yvect==newYVect[length(newYVect)])

	## strip off both unneeded X and Y values
	sliceMask <- sliceMask[minX:maxX, minY:maxY]

	## list contains the slice mask, with the new xvect and yvect
	## data spans (all coordinates existing in the slice)
	return (list(sliceMask=sliceMask, xIDx=minX:maxX, yIDx=minY:maxY))
}


## .getFunctionResults------------------2015-04-10
## function takes in a list of sliceValues which is a list of data vectors
## for each polygon in sliceValues that has not been clipped out due to xlim
## and or ylim arguments. This function applies every function to the list values
## and puts the result into a data frame for each polygon function combination.
## Function returns a list of data frames.
## ---------------------------------------------NL
.getFunctionResults <- function(sliceValues, functions, validPolys, na.rm)
{
	## run on every list element (polygon) of the sliceValues list
	functionResultsList <- lapply(sliceValues, function(sliceValue, fun, na.rm)
	{
		## apply all functions to each list element (polygon), returns list of function values.
		## functionResults is a list of all results from functions for an individual polygon
		## unlist is called to make functionResults a vector of returned function values.
		functionResults <- list()
		for(functionNumber in 1:length(functions)){
			if(length(na.rm) > 1) {
				## call all functions on a polygon and return the result
				if (na.rm[functionNumber]) {
					cmd <- paste("result <- ", functions[functionNumber], "(na.omit(sliceValue))", sep="")
					expression <- parse(text=cmd)
				} else {
					cmd <- paste("result <- ", functions[functionNumber], "(sliceValue)", sep="")
					expression <- parse(text=cmd)
				}
			} else {
				if(na.rm) {
					cmd <- paste("result <- ", functions[functionNumber], "(na.omit(sliceValue))", sep="")
					expression <- parse(text=cmd)
				}
				else {
					cmd <- paste("result <- ", functions[functionNumber], "(sliceValue)", sep="")
					expression <- parse(text=cmd)
				}
			}
			result <- eval(expression)
			functionResults <- c(functionResults, result)
		}
		## converts a list of function values to a vector of function values
		return (unlist(functionResults))
	}, functions, na.rm)

	## converts list of polygon, with vectors of function values, into a singular vector
	vectorOfResults <- unlist(functionResultsList)

	## convert singular vector into a matrix
	dim(vectorOfResults) <- c(length(functions), length(functionResultsList))
	matrixOfResults <- vectorOfResults

	## transpose matrix to create a proper data frame
	matrixOfResults <- t(matrixOfResults)

	## create dataframe from matrix
	df <- data.frame(validPolys, matrixOfResults)

	## gives function names to columns
	colnames(df) <- c("PID", functions)
	return(df)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.getFunctionResults


## .getNumberSequenceInfo---------------2015-04-10
## function takes in a locationSet (returned from findPolys function)
## and a combine value. Returns a list containing a numberSequence vector that
## is used for clipping vectors properly with a combine argument, and a
## validPolys vector that contains only the polygon ID's that are still in range
## after the user has specified a xlim/ylim clip.
## ---------------------------------------------NL
.getNumberSequenceInfo <- function(locationSet, combine)
{
	numberSequence <- NULL
	if (combine > 0) {
		numberSequence <- rep(sort(locationSet$PID), each=combine)
	}
	validPolys <- sort(unique(locationSet$PID))
	return (list(numberSequence=numberSequence, validPolys=validPolys))
}


## .getPolysVectorList------------------2015-04-10
## function takes in ncdfData slices and returns a list with each polygons
## data points that fall inside the polygon region.
## ---------------------------------------------NL
.getPolysVectorList <- function(slice, polygonList)
{
	polygonDataList <- lapply(polygonList, function(polygon, slice)
	{
		return(slice$data[polygon[["EID"]]])
	}, slice)
	return(polygonDataList)
}


.getSlice <- function(ncdfData, slice)
{
	if (!is.numeric(slice) && !is.character(slice)) {
		stop("Slice must be either an integer value or an existing date string in ncdfData")
	}
	## attempt to extract the slice
	if (is.numeric(slice)) {
		slice <- as.integer(slice)
		if (slice < 1 || slice > length(ncdfData)) {
			stop("Unable to extract slice: integer does not index an existing slice.")
		}
	}
	slice <- ncdfData[[slice]]

	## user did not specify a correct/existing slice; integer case handled above
	if (is.null(slice)) {
		stop("Unable to extract slice: date string not a valid slice name.")
	}
	return (slice)
}


## .getSliceNames-----------------------2018-10-23
## function converts a vector of time stamps 
## to YYYY-MM-DD HH:MM:SS format
## ------------------------------------------NL|RH
.getSliceNames <- function(times, unit, origin)
{
	if (grepl("^[Yy]ears$", unit)){
		times <- as.character(as.POSIXlt(times * (365.25 * 24 * 60 * 60), tz="GMT", origin=origin))
	} else if (grepl("^[Mm]onths$", unit)){
		times <- as.character(as.POSIXlt(times * (30.4375 * 24 * 60 * 60), tz="GMT", origin=origin))
	} else if (grepl("^[Ww]eeks$", unit)){
		times <- as.character(as.POSIXlt(times * (7 * 24 * 60 * 60), tz="GMT", origin=origin))
	} else if (grepl("^[Dd]ays$", unit)){
		times <- as.character(as.POSIXlt(times * (24 * 60 * 60), tz="GMT", origin=origin))
	} else if (grepl("^[Hh]ours$", unit)){
		times <- as.character(as.POSIXlt(times * (60 * 60), tz="GMT", origin=origin))
	} else if (grepl("^[Mm]inutes$", unit)){
		times <- as.character(as.POSIXlt(times * 60, tz="GMT", origin=origin))
	} else if(grepl("^[Ss]econds$", unit)){
		times <- as.character(as.POSIXlt(times, tz="GMT", origin=origin))
	} else {
		warning("Was not able to properly label slice names")
	}
	return (times)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.getSliceNames


.ncdfDataTlimClip <- function(ncdfData, tlim, nc=NULL)
{
	## save slice names
	sliceNames <- names(ncdfData)
	## remove slice names from ncdfData
	## when reassigning names old cause trouble with reassignment 
	names(ncdfData) <- NULL 
	
	## backup of ncdfData attributes
	attributesNcdfData <- attributes(ncdfData)
	if(!is.character(tlim)){
		stop("tlim must be of type character please specify a date string e.g., 2000-01-01")
	}
	if (length(tlim) != 2) {
		stop("tlim must be a vector of two strings, each formatted as \"YYYY-MM-DD HH:MM:SS\"")
	}
	## convert date strings into POSIXlt objects for date comparison
	dateNames <- as.POSIXlt(sliceNames, tz="GMT")

	## keep only the dates that fall in tlim
	datesInTLIM <- dateNames >= as.POSIXlt(tlim[1], tz="GMT") &
	dateNames <= as.POSIXlt(tlim[2], tz="GMT")

	## clip ncdfData to only keep slices in tlim
	ncdfData <- ncdfData[datesInTLIM]
	attributes(ncdfData) <- attributesNcdfData
	names(ncdfData) <- sliceNames[datesInTLIM]
	return(ncdfData)
}


## .ncdfDataXClip-----------------------2015-04-10
## function takes in a ncdfData object and removes data that 
## does not fall in the specified xlim region
## ---------------------------------------------NL
.ncdfDataXClip <- function(ncdfData, xlim, include.lowest)
{
	if(class(xlim) != "numeric"){
		stop("xlim must be a numeric range vector of length two")
	}
	if(length(xlim) != 2){
		stop("xlim must contain two numeric values")
	}
	if(xlim[1] > xlim[2]){
		stop("xlim must vectors must be increasing e.g., c(min,max)")
	}
	## get x vector attribute containing the longitude information from ncdfData
	xvect <- attr(ncdfData, "x")

	## ensure proper boundary placement with include.lowest argument
	if(is.null(include.lowest)) {
		newXVect <- xvect[xvect <= xlim[2] & xvect >=  xlim[1]]
	} else if(include.lowest==TRUE) {
		newXVect <- xvect[xvect < xlim[2] & xvect >=  xlim[1]]
	} else if(include.lowest==FALSE) {
		newXVect <- xvect[xvect <= xlim[2] & xvect >  xlim[1]]
	} else {
		stop("include.lowest must be of value TRUE, FALSE, or NULL")
	}
	minX <- which(xvect==newXVect[1])
	maxX <- which(xvect==newXVect[length(newXVect)])

	## get x vector attribute containing the longitude information from ncdfData
	attr(ncdfData, "x") <- newXVect

	## make a backup of attributes, apply them after lapply function
	ncdfDataAttributes <- attributes(ncdfData)

	ncdfData <- lapply(ncdfData, function(slice, minX, maxX)
	{
		for(layerName in names(slice)){
			slice[[layerName]] <- slice[[layerName]][minX:maxX,] 
		}
		return(slice)
	}, minX, maxX)
	attributes(ncdfData) <- ncdfDataAttributes
	return(ncdfData)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.ncdfDataXClip


## .ncdfDataYClip-----------------------2015-04-10
## function takes in a ncdfData object and removes data that
## does not fall in the specified ylim region
## ---------------------------------------------NL
.ncdfDataYClip <- function(ncdfData, ylim, include.lowest)
{
	if(class(ylim) != "numeric"){
		stop("ylim must be a numeric range vector of length two")
	}
	if(length(ylim) != 2){
		stop("ylim must contain two numeric values")
	}
	if(ylim[1] > ylim[2]){
		stop("ylim must vectors must be increasing e.g., c(min,max)")
	}
	## get y vector attribute containing the latitude information from ncdfData
	yvect <- attr(ncdfData, "y")

	## ensure proper boundary placement with include.lowest argument
	if(is.null(include.lowest)) {
		newYVect <- yvect[yvect <= ylim[2] & yvect >=  ylim[1]]
	} else if(include.lowest==TRUE) {
		newYVect <- yvect[yvect < ylim[2] & yvect >=  ylim[1]]
	} else if(include.lowest==FALSE) {
		newYVect <- yvect[yvect <= ylim[2] & yvect >  ylim[1]]
	} else {
		stop("include.lowest must be of value TRUE, FALSE, or NULL")
	}
	minY <- which(yvect==newYVect[1])
	maxY <- which(yvect==newYVect[length(newYVect)])

	attr(ncdfData, "y") <- newYVect

	## make a backup of attributes, apply them after lapply function
	ncdfDataAttributes <- attributes(ncdfData)

	ncdfData <- lapply(ncdfData, function(slice, minY, maxY)
	{
		for(layerName in names(slice)){
		  slice[[layerName]] <- slice[[layerName]][,minY:maxY] 
		}
		return(slice)
	}, minY, maxY)
	attributes(ncdfData) <- ncdfDataAttributes
	return(ncdfData)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.ncdfDataYClip


## .pInv--------------------------------2018-11-09
## Find nearest position in vector choice using a target point.
## source: ## https://stat.ethz.ch/pipermail/r-help/2008-July/167216.html
## ---------------------------------------------RH
.pInv = function(p,v){
	## Using sapply allows multiple target points p
	sapply(p, function(x,v){
		 ## occasionally two vector points are equidistant to the target p
		which(abs(v-x)==min(abs(v-x)))[1]
	}, v=v)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.pInv
