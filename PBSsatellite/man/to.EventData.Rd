\name{to.EventData}
\alias{to.EventData}
\title{Convert ncdfData Slice to EventData}
\description{
  Create \code{\link{EventData}} object from a \code{\link{ncdfData}} slice. 
}
\usage{
  to.EventData(ncdfData, slice)
}
\arguments{
  \item{ncdfData}{\code{\link{ncdfData}} where \code{slice} is located.  
    (required).}
  \item{slice}{date string or integer of \code{slice} location.}
}
\details{
  Converts a \code{\link{ncdfData}} \code{slice} to \code{\link{EventData}}.
  \code{\link{EventData}} makes \code{\link{ncdfData}}
  compatible with \code{\link{PBSmapping}} functionality. \code{\link{EventData}}
  is used to find which data points are in a polygon and which points
  fall outside a polygon, known as the points in polygon problem. 
}
\value{
  \code{\link{EventData}} with \code{\link{ncdfData}} \code{slice} information. 
}
\seealso{
  \code{\link{extractTimeSeries}}, \code{\link{assessMissingData}},
  \code{\link{clipRegion}}, \code{\link{PBSmapping}}, 
  \code{\link{findPolys}}.
}
\examples{
local(envir=.PBSsatEnv,expr={
  ## load ncdfData object
  data(sst)
  
  ## convert slice to ncdfData
  ed <- to.EventData(sst, slice=1)
  dim(ed); head(ed)
})
}
\author{
  \href{https://ca.linkedin.com/in/nicholas-lefebvre-b94b21132}{Nicholas Lefebvre}, Software Developer, Edmonton AB\cr
  Last modified \code{Rd: 2016-03-23}
}
\keyword{manip}
