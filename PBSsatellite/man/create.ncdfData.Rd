\name{create.ncdfData}
\alias{create.ncdfData}
\title{
  Create a NetCDF File
}
\description{
  Create a NetCDF (ncdf4) file from data.
}
\usage{
create.ncdfData(filename, xvals, yvals, tvals, 
   tmess="days since 1900-01-01", zlist, mv=-99, 
   dataType="Chlorophyll", dataUnits="mg/m3", longname=filename)
}
\arguments{
  \item{filename}{\code{character} -- name of NetCDF file to create.}
  \item{xvals}{\code{numeric} -- discrete values of X (longitude) that define the spatial grid.}
  \item{yvals}{\code{numeric} -- discrete values of Y (latitude) that define the spatial grid.}
  \item{tvals}{\code{character} -- discrete values of time (dates \code{"YYYY-MM-DD"}) to define temporal slices.}
  \item{tmess}{\code{character} -- time message or descriptor -- most versatile = \code{"days since YYYY-MM-DD"}, where specified date is the first in the series.}
  \item{zlist}{\code{numeric} -- list of z-value 3-D arrays, where each array dimension is defined by (lon, lat, date).}
  \item{mv}{\code{numeric} -- missing value indicator, usually -99 in \bold{PBSsatellite}.}
  \item{dataType}{\code{character} -- string representing type of data in the file (e.g., \code{"SST"}, \code{"Chl"}).}
  \item{dataUnits}{\code{character} -- string representing units of \code{dataType} (e.g., \code{"Celsius"}, \code{"mg m^3"}).}
  \item{longname}{\code{character} -- long name descriptor -- e.g., \dQuote{Some groovy chlorophyll data I downloaded from NOAA}.}
}
\details{
  This function creates NetCDF files in ncdf4 format, which loads into R quickly.
  
}
\value{
  An NetCDF file (e.g., \code{"some_chla.nc"}) is created in the user's working directory. 
  No object is returned to the user's working environment.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2018-10-24}
}
\note{
  See posting by by user3710546  (Mar 11 '15 at 10:48):\cr
  \href{https://stackoverflow.com/questions/28949971/writing-data-to-a-netcdf-file-with-r}{Writing data to a netCDF file with R}
}
\seealso{
  \code{\link[PBSsatellite]{convert.ncdfData}},
  \code{\link[PBSsatellite]{read.ncdfData}},
  \code{\link[PBSsatellite]{print.ncdfData}},
  \code{\link[PBSsatellite]{plot.ncdfData}}
}
\keyword{manip}

