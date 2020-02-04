\name{plot.ncdfData}
\alias{plot.ncdfData}
\title{Plot ncdfData Slice}
\description{
  Plot an \code{ncdfData} time slice.
}
\usage{
  \method{plot}{ncdfData}(x, slice, layer="data",
   xlim=NULL, ylim=NULL, style=c("image","contour"), projection="LL",
   tck=-0.014, tckMinor=0.5*tck, showXY=FALSE, ...)
}
\arguments{
  \item{x}{\code{\link{ncdfData}} -- object with multiple time slices.}
  \item{slice}{\code{numeric} -- time slice to plot; if \code{NULL} then first
    \code{slice} is selected}
  \item{layer}{\code{character} -- layer name to plot.}
  \item{xlim}{\code{numeric} -- range of X-coordinates to plot.}
  \item{ylim}{\code{numeric} -- range of Y-coordinates to plot.}
  \item{style}{\code{character} -- method for plotting the Z-value -- either \code{"image"} or \code{"contour"}.}
  \item{projection}{\code{character} -- desired projection when \code{\link{PolySet}} lacks
    a \code{projection} attribute; one of \code{"LL"}, \code{"UTM"}, or a
    numeric value. If Boolean, specifies whether to check \code{polys} for
    a \code{projection} attribute.}
  \item{tck}{\code{numeric} -- vector (length 1 or 2) describing the length of
    tick marks as a fraction of the smallest dimension. 
    If \code{tckLab==TRUE}, these tick marks will be automatically labelled.
    If given a two-element vector, the first element describes the tick marks on the
    x-axis and the second element describes those on the y-axis.}
  \item{tckMinor}{\code{numeric} -- vector (length 1 or 2) describing the length
    of tick marks as a fraction of the smallest dimension. These tick
    marks can not be automatically labelled. If given a two-element
    vector, the first element describes the tick marks on the x-axis and
    the second element describes those on the y-axis.}
  \item{showXY}{\code{logical} -- if \code{TRUE}, show points that represent each grid cell.}
  \item{...}{additional arguments sent to \code{style} function.}
}
\details{
  Plots an \code{\link{ncdfData}} \code{layer}. If no \code{slice} is given,
  the code selects the first \code{slice}. If no layer is given, the code uses
  \code{"data"} layer of \code{\link{ncdfData}} object. 
  
  The user can select a region to plot based on \code{xlim} and/or 
  \code{ylim} arguments.
  
  The user can select different plotting styles to plot \code{\link{ncdfData}}
  such as \code{"image"}, or \code{"contour"}.
}
\seealso{
  \code{\link[PBSsatellite]{ncdfData}},
  \code{\link[PBSsatellite]{create.ncdfData}},
  \code{\link[PBSsatellite]{read.ncdfData}}\cr
  \bold{PBSmapping:}\cr
  \code{\link[PBSmapping]{plotMap}},
  \code{\link[PBSmapping]{addPolys}}.
}
\note{
  If an \code{ncdfData} object is created from irregular (X,Y) coordinates, 
  grid cells are also irregular with boundaries half the distance between 
  X-neighbours and between Y-neighbours.
}
\examples{
local(envir=.PBSsatEnv,expr={
  ## load ncdfData object
  data(sst)
  ## load worldLL polygons from PBSmapping
  data(worldLL)
  
  ## plot map using image functionality on the first slice
  plot(sst, slice=1, style="image")
  addPolys(worldLL)
  
  ## plot map using contour functionality on the first slice
  plot(sst, slice=1, style="contour")
  addPolys(worldLL)
})
}
\author{
  \href{https://academic.macewan.ca/boersn/}{Nicholas M. Boers}, Associate Professor -- Computer Science\cr
  MacEwan University, Edmonton AB

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2018-10-22}
}
\keyword{hplot}
