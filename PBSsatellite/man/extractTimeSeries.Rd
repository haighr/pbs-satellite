\name{extractTimeSeries}
\alias{extractTimeSeries}
\title{Extract ncdfData Time Series}
\description{
  Create a time series of \code{\link{ncdfData}} based on a given location, 
  using specified functions.
}
\usage{
extractTimeSeries (ncdfData, xlim=NULL, ylim=NULL,
                              polygons=NULL, functions=c("sum", "mean", "sd"),
                              na.rm=TRUE, tlim=NULL, combine=1, by=NULL, 
                              include.lowest=TRUE)
}
\arguments{
  \item{ncdfData}{\code{\link{ncdfData}} to extract time series from 
    (required).}
  \item{xlim}{range of X-coordinates.}
  \item{ylim}{range of Y-coordinates.}
  \item{polygons}{complex range of coordinates from which to extract time
    series. If more than one polygon time series is created for each polygon
    and will be specified by a PID.}
  \item{functions}{vector of function strings, where each function accepts
    a numeric vector and produce a single numeric value.}
  \item{na.rm}{Boolean or Boolean vector that indicates whether NA's
    should be omitted;
    if a vector, should match one-to-one with the functions.}
  \item{tlim}{start date and end date of when to begin and end a time series.}
  \item{combine}{integer number of slices to be combined into each
    call to each function in \code{functions}, e.g., given
    \code{\link{ncdfData}} with 6 slices, a \code{combine} value 
    of 2 will produce time series statistics for three times, where each of its      
    times considers the data from two source slices.}
  \item{by}{integer number indicating whether the function should skip slices,
    e.g., produce time series statistics for every second slice in
    \code{\link{ncdfData}}.}
  \item{include.lowest}{see \code{\link{clipRegion}}.}
}
\details{
  In the case of \code{xlim}/\code{ylim}, the resulting time
  series data does not include a subregion identifier.  In the case of
  \code{polygons} (one or more), the resulting data set contains a
  subregion identifier equal to the corresponding PID from
  \code{polygons}.
  
  For each slice in the data set, the function will determine which points
  fall within the region(s) of interest.  It will pass these points (as a
  vector) into each of the listed \code{functions}.
  
  If a \code{combine} value is provided that is not a factor of
  \code{length(ncdfData)} (# of slices), slices will be removed
  from the tail of the \code{\link{ncdfData}} object in order to
  accommodate the \code{combine} value.
}
\value{
  A list of date named data frames. The data frames contain a PID for
  each \code{polygon} that exists inside of \code{xlim} and
  \code{ylim}. If \code{xlim} and/or \code{ylim=NULL functions} will
  take the full span of \code{\link{ncdfData}}. Each polygon contains a
  row in the data frame and will have columns for every function in
  \code{functions}.

  If \code{polygons=NULL} there will be only one PID for the entire
  \code{xlim/ylim} of \code{\link{ncdfData}}.
}
\seealso{
  \code{\link{clipRegion}}, \code{\link{EventData}}.
}
\examples{
local(envir=.PBSsatEnv,expr={
  ## load ncdfData object
  data(sst)
  
  ## create a time series based on full map
  ts <- extractTimeSeries(sst)
  print(ts)
  
  ## create a time series based on 2 polygons
  ## create polygons
  polys <- data.frame(PID=c(rep(1, 4), rep(2, 4)), POS=c(1:4, 1:4),
                     X=c(155, 160, 150, 180, 0, 20, 20, 0),
                     Y=c(75, 50, 10, 85, 20, 20, 40, 40))
  ts <- extractTimeSeries(sst, polygons=polys)
  print(ts)
})
}
\author{
  \href{https://ca.linkedin.com/in/nicholas-lefebvre-b94b21132}{Nicholas Lefebvre}, Software Developer, Edmonton AB\cr
  Last modified \code{Rd: 2016-03-23}
}
\keyword{math}
