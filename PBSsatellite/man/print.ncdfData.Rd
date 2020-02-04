\name{print.ncdfData}
\alias{print.ncdfData}
\title{Print ncdfData Object}
\description{
  Pring ncdfData object. 
}
\usage{
 \method{print}{ncdfData}(x, \dots)
}
\arguments{
  \item{x}{\code{\link{ncdfData}} object to be printed (required).}
  \item{...}{additional printing arguments.}
}
\details{
  Prints a \code{\link{ncdfData}} object.
}
\seealso{
  \code{\link{read.ncdfData}}, \code{\link{ncdfData}}, \code{\link{print}}.
}
\examples{
local(envir=.PBSsatEnv,expr={
  ## load ncdfData object
  data(sst)
  
  ## print ncdfData object
  print(sst)
})
}
\author{
  \href{https://academic.macewan.ca/boersn/}{Nicholas M. Boers}, Associate Professor -- Computer Science\cr
  MacEwan University, Edmonton AB\cr
  Last modified \code{Rd: 2016-03-23}
}
\keyword{print}
