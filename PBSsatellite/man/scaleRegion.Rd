\name{scaleRegion}
\alias{scaleRegion}
\title{Scale ncdfData to a New Resolution}
\description{
  Scale ncdfData slices to a new resolution based on a scale factor.
}
\usage{
scaleRegion(ncdfData, scaleFactor, fun="drop", placement="topleft",
      includeErrorMatrix=FALSE, includeMissMatrix=FALSE, 
      remainder="crop", na.rm=TRUE)
}
\arguments{
  \item{ncdfData}{\code{\link{ncdfData}} which will be scaled by
    \code{scaleFactor} (required).}
  \item{scaleFactor}{positive or negative integer describing the scale factor.
    A positive integer will scale up a \code{\link{ncdfData}} object a negative
    integer will scale down a \code{\link{ncdfData}} object. All integers must 
    be a power of two. A positive integer increases the number of data points 
    to 1*\code{scaleFactor} in each axis for a total increase of
    1*\code{scaleFactor}^2. A negative integer reduces the number
    of data points to 1/\code{scaleFactor} in each axis for a
    total reduction of 1
    /\code{scaleFactor}^2 (required).}
  \item{fun}{string of a function to used to scale down: \code{"mean"},
    \code{"min"}, \code{"max"}, \code{"drop"}. When
    scaling up \code{"repeat"} is always used.}
  \item{placement}{string indicating placement for the computed data point:
    \code{"topleft"}, \code{"centre"}.}
  \item{includeErrorMatrix}{logical indicates whether an error matrix should
    paired with each data matrix in the resulting \code{\link{ncdfData}}
    object.}
  \item{includeMissMatrix}{logical indicates whether a missing matrix should be
    paired with each data matrix in the resulting
    \code{\link{ncdfData}} object.}
  \item{remainder}{string if \code{"crop"}: if \code{len(x)} of
    \code{\link{ncdfData}} and/or \code{len(y)} of \code{\link{ncdfData}}
    is not a factor of \code{scaleFactor} crop will remove rows and/or columns
    from \code{\link{ncdfData}} slices in order to make \code{len(x)}
    and/or \code{len(y)} a factor of \code{\link{ncdfData}}. If
    \code{"fill"}: rows and/or columns of NA values will be added on to
    \code{\link{ncdfData}} slices to make \code{len(x)} and/or
    \code{len(y)} a factor of \code{scaleFactor}.}
  \item{na.rm}{logical indicates if vector values should omit NA values before 
  call to \code{fun}.}          
}
\details{
  It is common for satellite data to be in different resolutions. 
  e.g., SST (1/4 degree) or Chl (1/8 degree)
  It is much easier to compare different data sets that are in a standardized 
  resolution. 
  
  This function creates a new \code{\link{ncdfData}} object with slices
  converted to a new resolution. For this initial version, the new
  resolution must be an integer scaleFactor of the original resolution.
  The user must specify a \code{scaleFactor}. When computing new
  data points, the user may choose to have the computed data point placed
  at the \code{"topleft"} point's position or in the
  \code{"centre"} of the scaled points (scaling down only). 
  
  A negative \code{scaleFactor} argument can take a \code{fun}
  to perform the scaling operation. Drop is also an option that
  \code{\link{scaleRegion}} provides. \code{fun="drop"}: drop
  points that do not fall on the points of the new scaled down region.
  
  A positive \code{scaleFactor} argument will only use the
  \code{"repeat"} method for version one.
  \code{fun="repeat"}: is the default and only function available
  for scaling up, repeat will repeat a points data 1 *
  \code{scaleFactor}^2 times in order to properly increase the
  scale of \code{\link{ncdfData}} slices.
}
\value{
  A \code{\link{ncdfData}} object containing slices with a newly scaled
  region. 
  If the user has specified \code{includeErrorMatrix=TRUE} and/or
  \code{includeMissMatrix=TRUE} the slices in the
  \code{\link{ncdfData}} object will now have an additional two layers
  (error and/or missing). These layers are of the same resolution (point
  for point) as the slices. These additional layers will have a percent
  error and a percentage of missing values with every corresponding
  point in a given slice. 
  
  If \code{placement="topleft"} \code{x} and \code{y}
  attributes will be increased (+ \code{scaleFactor}) or reduced
  (- \code{scaleFactor}).
  
  If \code{placement="centre"} \code{x} and \code{y} 
  attributes will be increased (+ \code{scaleFactor}) or reduced
  (- \code{scaleFactor}) and \code{x} and \code{y}
  attributes will be shifted towards the centre of the point.
}

\seealso{
  \code{\link{read.ncdfData}}.
}
\examples{
local(envir=.PBSsatEnv,expr={
  ## load ncdfData object
  data(sst)
  
  ## scale down ncdfData slices by a factor of 2, using mean
  sd <- scaleRegion(sst, scaleFactor=-2, fun="mean", remainder="fill")
  print(sd)
  
  ## scale down ncdfData slices by a factor of 2, using drop, place result
  ## in the centre of the clip region
  sd2 <- scaleRegion(sst, scaleFactor=-2, fun="drop", remainder="drop",
  placement="centre")
  print(sd2)
  
  ## scaling up ncdfData slices by a factor of 4
  sd3 <- scaleRegion(sst, scaleFactor=4, fun="repeat")
  print(sd3)
})
}
\author{
  \href{https://ca.linkedin.com/in/nicholas-lefebvre-b94b21132}{Nicholas Lefebvre}, Software Developer, Edmonton AB\cr
  Last modified \code{Rd: 2016-03-23}
}
\keyword{manip}
