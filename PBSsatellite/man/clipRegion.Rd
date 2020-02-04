\name{clipRegion}
\alias{clipRegion}
\title{Clip an existing ncdfData object}
\description{
 Clip the region of a \code{\link{ncdfData}} object. Clipping is applied to
 all slices in \code{\link{ncdfData}} for the specified region.
}
\usage{
clipRegion(ncdfData, xlim=NULL, ylim=NULL, polygons=NULL, include.lowest=TRUE) 
}
\arguments{
  \item{ncdfData}{\code{\link{ncdfData}} from which to clip region (required).}
  \item{xlim}{range of X-coordinates for data slice(s).}
  \item{ylim}{range of Y-coordinates for data slice(s).}
  \item{polygons}{data slice(s) dimensions will be created to match one or more
    \code{polygon} dimensions.}
  \item{include.lowest}{ignored unless user specifies
    \code{xlim/ylim}: if \code{TRUE}, includes points that
    fall on \code{min(xlim)} and/or \code{min(ylim)} but
    not points that fall on \code{max(xlim)} and/or
    \code{max(ylim)}. If \code{FALSE}, it does the
    opposite, includes points that fall on \code{max(xlim)}
    and/or \code{max(ylim)} but not points that fall on
    \code{min(xlim)} and/or \code{min(ylim)}.  If
    \code{NULL}, includes all boundary points.}
}
\details{
  In most cases data sets contain information spanning the whole
  world; therefore, it is useful for the user to be able select
  a region that better suits their needs.  The user will be able to
  select a geographical area by specifying \code{xlim} and/or
  \code{ylim} arguments, e.g., Northeast Pacific.  For more
  complex selection a user is able to select a geographical area based
  on a \code{polygon}, e.g., Georgia Strait.

  Clipping will be applied to every \code{\link{ncdfData}} slice.
  
  In situations where \code{xlim} and/or \code{ylim} as
  well as a polygon are provided, \code{xlim} and
  \code{ylim} clipping with be clipped first from
  \code{\link{ncdfData}}, the \code{polygon} clipping will be
  applied to the result. 
}
\value{
  The \code{\link{clipRegion}} method returns a new \code{\link{ncdfData}} 
  object containing geographically modified slices from an existing
  \code{\link{ncdfData}}.
}
\seealso{
  \code{\link{extractSlices}}.
}
\examples{
local(envir=.PBSsatEnv,expr={
  ## load ncdfData object
  data(sst)
  ## load worldLL polygons for displaying
  data(worldLL)
  
  ## clip region based on xlim and ylim
  ncdfDataClip <-clipRegion(sst, xlim=c(190, 320), ylim=c(5, 80), 
                                                   include.lowest=NULL)
  ## print newly clipped ncdfData object
  print(ncdfDataClip)
  
  ## plot ncdfData object
  plot(ncdfDataClip, slice=1)
  addPolys(worldLL, col="beige")
  
  ## clip region based on xlim, ylim, and polygons
  ## create 2 polygons
  polys <- data.frame(PID=c(rep(1, 4), rep(2, 4)), POS=c(1:4, 1:4),
                     X=c(155, 160, 150, 180, 0, 20, 20, 0),
                     Y=c(75, 50, 10, 85, 20, 20, 40, 40))
                     
  ncdfDataClip <-clipRegion(sst, xlim=c(.5, 300), ylim=c(-50, 90), 
  polygons=polys, include.lowest=NULL)
  
  ## print newly clipped ncdfData object
  print(ncdfDataClip)
  
  ## plot ncdfData object
  plot(ncdfDataClip, slice=1)

  ## add some polygons to show the clipped region
  addPolys(polys, border="blue", lwd=2)
  addPolys(worldLL, border="gray")
})
}
\author{
  \href{https://ca.linkedin.com/in/nicholas-lefebvre-b94b21132}{Nicholas Lefebvre}, Software Developer, Edmonton AB\cr
  Last modified \code{Rd: 2016-03-23}
}
\keyword{manip}
