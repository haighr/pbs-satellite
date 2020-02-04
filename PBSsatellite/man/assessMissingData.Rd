\name{assessMissingData}
\alias{assessMissingData}
\title{Assessment of Missing Data}
\description{
  Create an assessment of missing data in each 
  \code{\link{ncdfData}} slice. 
}
\usage{
assessMissingData(ncdfData, tlim=NULL, xlim=NULL, ylim=NULL, 
      polygons=NULL, include.lowest=TRUE)
}
\arguments{
  \item{ncdfData}{\code{\link{ncdfData}} used for missing data assessment 
    (required).}
  \item{tlim}{start date and end date of assessment. If
    \code{tlim=NULL}, then assess all slices.}
  \item{xlim}{range of X-coordinates.}
  \item{ylim}{range of Y-coordinates.}
  \item{polygons}{complex region to apply assessment.}
   \item{include.lowest}{see \code{\link{clipRegion}}. }
}
\details{
  It is common with satellite data for data sets to have missing
  values. This function is useful to indicate to the user whether a given
  data set is relatively complete or incomplete. It is also common for
  satellites to have a bad reading for a given duration, e.g., a week,
  thus the function is applied to every slice. The user can determine
  which slices are complete enough for their use.
}
\value{
  Numeric vector containing a percentage of missing data for each slice
  in \code{\link{ncdfData}}.
}
\seealso{
  \code{\link{extractSlices}}, \code{\link{clipRegion}}.
}
\examples{
local(envir=.PBSsatEnv,expr={
  ## load ncdfData object
  data(sst)
  ## get missing data assessment
  md <- assessMissingData(sst)
  print(md)
  
  ## user a polygon for missing data assessment
  ## create 2 polygons
  polys <- data.frame(PID=c(rep(1, 4), rep(2, 4)), POS=c(1:4, 1:4),
                     X=c(155, 160, 150, 180, 0, 20, 20, 0),
                     Y=c(75, 50, 10, 85, 20, 20, 40, 40))
  md <- assessMissingData(sst, polygons=polys)
  print(md)
})
}
\author{
  \href{https://ca.linkedin.com/in/nicholas-lefebvre-b94b21132}{Nicholas Lefebvre}, Software Developer, Edmonton AB\cr
  Last modified \code{Rd: 2016-03-23}
}
\keyword{math}
