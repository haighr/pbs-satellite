\name{convert.ncdfData}
\alias{convert.ncdfData}
\title{
  Convert ASCII Data to NetCDF
}
\description{
  Convert satellite ASCII data to a netCDF file, specifically a ncdf4 binary file.
}
\usage{
convert.ncdfData(filename, zfld, nc.filename="converted.nc",
   summary.func=sum, offset=c(0,0), mv=-99,
   dataType="Chlorophyll", dataUnits="mg/m3")
}
\arguments{
  \item{filename}{\code{character} -- name of the ASCII source file (comma-delimied, CSV-like) containing satellite gridded data.}
  \item{zfld}{\code{character} -- vector of field names in the source file that contain satellite measurements (e.g., \code{"Chl"}).}
  \item{nc.filename}{\code{character} -- name of the NetCDF binary file that the user wants to create.}
  \item{summary.func}{\code{function} -- summary function (e.g., \code{sum}) to aggregate multiple measurements at unique combinations of (lon, lat, time).}
  \item{offset}{\code{numeric} -- coordinate (x,y) offset to adjust ASCII (X,Y) grid coordinates in case the latter is defined by some vertex other than the top left one.}
  \item{mv}{\code{numeric} -- missing value indicator, usually -99 in \bold{PBSsatellite}.}
  \item{dataType}{\code{character} -- string representing type of data in the file (e.g., \code{"SST"}, \code{"Chl"}).}
  \item{dataUnits}{\code{character} -- string representing units of \code{dataType} (e.g., \code{"Celsius"}, \code{"mg m^3"}).}
}
\details{
  Users sometimes prefer storing satellite data in cumbersome ASCII files.
  This function attempts to convert such files to a more effienct ncdf4 binary format.
  The function imports the ASCII file and locates the appropriate three dimesnions (lon, lat, date) to create 3-D arrays of z-value data.
  These data are then passed to the \bold{PBSsatellite} function \code{'create.ncdfData'}, which creates an ncdf4 binary file.
}
\value{
  A NetCDF (ncdf4) file, e.g. \code{"converted_chla.nc"}, is created in the user's working directory. 
  No object is returned to the user's working environment.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2018-10-24}
}
\note{
  This function uses the \code{'fread'} function from the R package \bold{data.table} to facilitate loading very large ASCII files.
}
\seealso{
  \code{\link[PBSsatellite]{create.ncdfData}},
  \code{\link[PBSsatellite]{read.ncdfData}},
  \code{\link[PBSsatellite]{plot.ncdfData}}
}
\keyword{manip}

