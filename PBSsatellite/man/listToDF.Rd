\name{listToDF}
\alias{listToDF}
\title{Convert a List to a Data Frame}
\description{
  Convert a Lists such as a TimeSeries returned from
  \code{\link{extractTimeSeries}} into a singular
  data frame.
}
\usage{
listToDF(lst, newColumnName="names")
}
\arguments{
  \item{lst}{\code{\link{list}} to be converted to a
    data frame.}
  \item{newColumnName}{name of the new column in the
    data frame.}
}
\details{
  Converts a list to a data frame. A list contains
  a name for each list element. The list names are extracted
  from the list and are stored in a column named by the
  \code{newColumnName} argument. This simplifies a list by having
  both the list names and the list elements data component in a
  single row contained in a data frame.
}
\value{
  A data frame containing a lists names as part
  of the data frame located in the column named
  \code{newColumnName}.
}
\seealso{
  \code{\link{data.frame}}, \code{\link{extractTimeSeries}}.
}
\examples{
local(envir=.PBSsatEnv,expr={
  ## load ncdfData object
  data(sst)

  ## create a time series based on 2 polygons
  ## create polygons
  polys <- data.frame(PID=c(rep(1, 4), rep(2, 4)), POS=c(1:4, 1:4),
                    X=c(155, 160, 150, 180, 0, 20, 20, 0),
                    Y=c(75, 50, 10, 85, 20, 20, 40, 40))
  ts <- extractTimeSeries(sst, polygons=polys)

  ## turn TimeSeries list into a DataFrame
  tsDF <- listToDF(ts)
  print(tsDF)
})
}
\author{
  \href{https://academic.macewan.ca/boersn/}{Nicholas M. Boers}, Associate Professor -- Computer Science\cr
  MacEwan University, Edmonton AB\cr
  Last modified \code{Rd: 2016-03-23}
}
\keyword{manip}
