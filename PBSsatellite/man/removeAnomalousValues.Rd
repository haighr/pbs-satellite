\name{removeAnomalousValues}
\alias{removeAnomalousValues}
\title{Remove Anomalous Values from a ncdfData Object}
\description{
  Remove specified anomalies from every \code{\link{ncdfData}} slice.
}
\usage{
removeAnomalousValues(ncdfData, zlim)
}
\arguments{
  \item{ncdfData}{\code{\link{ncdfData}} from which to remove anomalies
    (required).}
  \item{zlim}{numeric vector containing a range of acceptable values in
    \code{\link{ncdfData}} slices. All values that do not fall in the
    range of \code{zlim} are removed from the data set and will be
    replaced with \code{\link{ncdfData}} objects missing value (required).
    NA can be used to omit part of a range in zlim.}
}
\details{
  It is common with satellite data for data sets to contain values that
  are anomalous. This can happen due to a variety of environmental
  reasons. A \code{zlim} argument is required that contains a
  numeric vector containing a range of values that are considered valid.
}
\value{
  returns a new \code{\link{ncdfData}} object containing slices with
  removed anomalous values.
}
\seealso{
  \code{\link{assessMissingData}}.
}
\examples{
local(envir=.PBSsatEnv,expr={
  ## load ncdfData object
  data(sst)
  
  ## remove values less than -2 and greater than 25
  newNcdfData <- removeAnomalousValues(sst, zlim=c(-2, 25))
  
  ## remove values greater than 25
  newNcdfData <- removeAnomalousValues(sst, zlim=c(NA, 25))
})
}
\author{
  \href{https://ca.linkedin.com/in/nicholas-lefebvre-b94b21132}{Nicholas Lefebvre}, Software Developer, Edmonton AB\cr
  Last modified \code{Rd: 2016-03-23}
}
\keyword{manip}
