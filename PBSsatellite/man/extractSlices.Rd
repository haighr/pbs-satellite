\name{extractSlices}
\alias{extractSlices}
\title{Extract ncdfData slice(s)}
\description{
 Extract slices from a \code{\link{ncdfData}} object.
}
\usage{
extractSlices(ncdfData, slices=NULL, dates=NULL, tlim=NULL)
}
\arguments{
  \item{ncdfData}{\code{\link{ncdfData}} from which to extract slices
    (required).}
  \item{slices}{numeric vector containing the slices to extract from
    \code{\link{ncdfData}}.}
  \item{dates}{vector of date strings containing specific dates of slices used 
    to extract slices from \code{\link{ncdfData}}.}
  \item{tlim}{vector of date strings (length 2) used to extract slices from 
  \code{\link{ncdfData}} based on a time range.}
}
\details{
  User must specify one of the three arguments \code{slices},
  \code{dates}, or \code{tlim}. A new \code{\link{ncdfData}}
  object will be created containing the slices that the user has indicated 
  using one of the three extraction methods. 
}
\value{
  A \code{\link{ncdfData}} object containing slices provided by the
  extraction method. 
}
\seealso{
  \code{\link{assessMissingData}}, \code{\link{clipRegion}}.
}
\examples{
local(envir=.PBSsatEnv,expr={
  ## load ncdfData object
  data(sst)
  
  ## extract slices based on date strings
  ## create dates object
  dates <- c("2001-03-01", "2001-05-01")
  newNcdfData <- extractSlices(sst, dates=dates)
  print(newNcdfData)
  
  ## extract slices based date range/tlim
  ## create tlim object
  tlim <- c("2001-02-04", "2001-07-02")
  newNcdfData <- extractSlices(sst, tlim=tlim)
  print(newNcdfData)
  
  ## extract slices based on slices
  newNcdfData <- extractSlices(sst, slices=c(2,3))
  print(newNcdfData)
})
}
\author{
  \href{https://ca.linkedin.com/in/nicholas-lefebvre-b94b21132}{Nicholas Lefebvre}, Software Developer, Edmonton AB\cr
  Last modified \code{Rd: 2016-03-23}
}
\keyword{manip}
