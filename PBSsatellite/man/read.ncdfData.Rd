\name{read.ncdfData}
\alias{read.ncdfData}
\alias{merge.ncdfData}
\title{
  Create New or Merged ncdfData Object
}
\description{
  Create and return an \code{ncdfData} object. When possible, auto
  detect names from file and inform the user. User has the ability to
  override any inconsistencies between NetCDF attribute names and the
  given \code{ncdfData} object attribute values.
}
\usage{
read.ncdfData(filename, 
   dataVariable=1, convertMissingValues=FALSE, 
   dataType=NULL, dataUnits=NULL, dataVname=NULL,
   xlim=NULL, ylim=NULL, tlim=NULL,
   x=NULL, y=NULL, time=NULL,
   Ux=NULL, Uy=NULL, Utime=NULL, verbose=FALSE) 

merge.ncdfData <- function (filenames,
   dataVariable=1, convertMissingValues=FALSE,
   dataType=NULL, dataUnits=NULL, dataVname=NULL,
   xlim=NULL, ylim=NULL, tlim=NULL,
   x=NULL, y=NULL, time=NULL,
   Ux=NULL, Uy=NULL, Utime=NULL, verbose=FALSE,
   read.ncdf=TRUE, write.ncdf=FALSE, outname)
}
\arguments{
  \item{filename, filenames}{\code{character} -- path/filename of the NetCDF file(s) to be read (required);
    \code{filenames} can also specific a collection of \code{ncdfData} object names when using \code{merge.ncdfData}.}
  \item{dataVariable}{\code{numeric} -- location of data variable within the NetCDF file \code{(default=1}).}
  \item{convertMissingValues}{\code{logical} -- if \code{TRUE}, convert missing
    values \code{"NA"} back to native form (\code{default=FALSE}).}
  \item{dataType}{\code{character} -- string representing type of data in the file (e.g.,\code{"SST"}, \code{"Chl"}).}
  \item{dataUnits}{\code{character} -- string representing units of dataType (e.g., \code{"Celsius"}, \code{"mg m^3"}).}
  \item{dataVname}{\code{character} -- name of the variable selected by \code{dataVariable} position;
    set \code{verbose=TRUE} to see what variables exist in NetCDF file(s) being read.}
  \item{tlim}{\code{character} -- range of time (slices) to import (use date strings like \code{"2000-01-01"});
    if \code{tlim=NULL} then import all slices.}
  \item{xlim}{\code{numeric} -- range of X-coordinates for data slice(s).}
  \item{ylim}{\code{numeric} -- range of Y-coordinates for data slice(s).}
  \item{x}{\code{character} -- name of x variable. If \code{xlim=NULL}, then \code{x="lon"}.}
  \item{y}{\code{character} -- name of y variable. If \code{ylim=NULL}, then \code{y="lat"}.}
  \item{time}{\code{character} -- name of time variable. If \code{time=NULL} then \code{time="time"}.}
  \item{Ux}{\code{character} -- units of x variable.}
  \item{Uy}{\code{character} -- units of y variable.}
  \item{Utime}{\code{character} -- units of time variable.}
  \item{verbose}{\code{logical} -- if \code{TRUE}, report the variables and dimensions within the NetCDF file being read.}
  \item{read.ncdf}{\code{logical} -- if \code{TRUE}, read in the specified NetCDF files to be merged;
    if \code{FALSE}, use existing \code{ncdfData} objects in user's working environment.}
  \item{write.ncdf}{\code{logical} -- if \code{TRUE}, write the merged \code{ncdfData} object to a NetCDF file \code{'.nc'}.}
  \item{outname}{\code{character} -- name of NetCDF output file \code{'.nc'}.}
}
\details{
  Creates an \code{ncdfData} object that can be used with other
  \bold{PBSsatellite} functions. Where possible, this function attempts to
  read names from the data file, but it allows the user to override
  names to account for the inconsistencies between different NetCDF files.
  The function \code{read.ncdfData} makes it possible for a
  variety of different NetCDF formats with varying data types to become
  compatible.

  Use \code{merge.ncdfData} to merge a series of similarly dimensioned
  NetCDF files or \code{ncdfData} objects. At present, the merge does not 
  check whether the dimensions among objects is consistent but assumes that the
  user has chosen files or objects from a series\cr e.g., 
  \href{https://www.ncdc.noaa.gov/data-access/marineocean-data/extended-reconstructed-sea-surface-temperature-ersst-v4}{Extended Reconstructed Sea Surface Temperature (ERSST) v4}.
}
\value{
  An \code{ncdfData} object containing attributes and data slices from 
  a one or more NetCDF files.
}
\seealso{
  \code{\link[PBSsatellite]{ncdfData}},
  \code{\link[PBSsatellite]{create.ncdfData}},
  \code{\link[PBSsatellite]{plot.ncdfData}},
  \code{\link[PBSsatellite]{clipRegion}},
  \code{\link[PBSsatellite]{extractSlices}}
}
\examples{
local(envir=.PBSsatEnv,expr={
  ## read in the whole NetCDF file containing the full region
  path <- system.file("sst.ltm.1971-2000.nc", package="PBSsatellite")
  ncdfData <- read.ncdfData(filename=path, verbose=TRUE)
  print(ncdfData)

  ## clipping the NetCDF file by dates
  ## create a tlim argument of date strings
  dates <- c("1-02-01","1-05-01")
  ncdfData <- read.ncdfData(filename=path, tlim=dates)
  print(ncdfData)

  ## clipping the NetCDF file by dates and region
  ncdfData <- read.ncdfData(filename=path, tlim=dates, xlim=c(20, 80),
    ylim=c(-80, 10))
  print(ncdfData)
})
}
\author{
  \href{https://ca.linkedin.com/in/nicholas-lefebvre-b94b21132}{Nicholas Lefebvre}, Software Developer, Edmonton AB\cr

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2018-10-24}
}
\keyword{manip}
