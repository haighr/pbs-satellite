\name{ncdfData}
\alias{ncdfData}
\title{ncdfData Object}
\description{
  PBS Satellite data object that contains satellite data for varying data types
  and spatial resolutions.
}

\details{
  An \code{\link{ncdfData}} object contains at least once slice. A slice
  contains a date name and at least one layer known as the data
  layer. The data layer contains the data information for gridded
  satellite data in a matrix format. The \code{\link{ncdfData}} objects
  contain attributes for data type (title of data), vectors of x and y
  coordinates, slice names, and data units.  Slices can optionally hold
  additional layers of information that contain point for point the same
  data span as the the data layer. Additional layers are created with
  functions such as \code{\link{scaleRegion}} that gives the user the
  option to include a missing layer and error layer when scaling down.
  
  A clipped \code{\link{ncdfData}} object with a polygon creates a
  complex region. In order to store a matrix style of representation of
  this complex region, slices layers use NaN as a place holder for data
  that has been clipped for the existing \code{\link{ncdfData}}
  object. Data that is missing from the original data set and has not
  been clipped is represented as NA.
  
}
\value{
  An \code{\link{ncdfData}} object
}
\seealso{
  \code{\link{read.ncdfData}}, \code{\link{scaleRegion}},
  \code{\link{clipRegion}}.
}
\author{
  \href{https://ca.linkedin.com/in/nicholas-lefebvre-b94b21132}{Nicholas Lefebvre}, Software Developer, Edmonton AB\cr
  \href{https://academic.macewan.ca/boersn/}{Nicholas M. Boers}, Associate Professor -- Computer Science\cr
  MacEwan University, Edmonton AB\cr
  Last modified \code{Rd: 2015-04-13}
}
\keyword{classes}
