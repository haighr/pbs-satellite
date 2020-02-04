\name{write.ncdfData}
\alias{write.ncdfData}
\title{
  Write ncdfData Slices to CSV Files
}
\description{
  Write specified slices from an \code{'ncdfData'} class object to comma-delimited text files.
  Each slice specified is written to a separate file.
}
\usage{
write.ncdfData(ncdfData, layer="data", slice, 
  file.prefix, save.dir="./ncdfData")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ncdfData}{\code{list} -- an object with class \code{'ncdfData'}.}
  \item{layer}{\code{character|numeric} -- layer in ncdfData object, usually \code{"data"}.}
  \item{slice}{\code{numeric|character} -- slice number(s) or name(s) to extract and dump to text file(s).}
  \item{file.prefix}{\code{character} -- prefix for output file names; if missing, 
    the code uses the \code{ncdfData}'s attribute called \code{'dataVname'}.}
  \item{save.dir}{\code{character} -- directory in which to save the slices of data as text files.}
}
\details{
  A simple function that recursively goes through the user-specified slices of an \code{ncdfData}
  object, creates matrices of Z values named by the attributes \code{'x'} and \code{'y'}, and
  writes them to comma-separated (\code{'.csv'}) files.
}
\value{
  Returns nothing, just dumps text files.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2018-11-09}
}
\seealso{
  \code{\link[PBSsatellite]{convert.ncdfData}},
  \code{\link[PBSsatellite]{create.ncdfData}},
  \code{\link[PBSsatellite]{download.ncdfData}},
  \code{\link[PBSsatellite]{merge.ncdfData}},
  \code{\link[PBSsatellite]{plot.ncdfData}},
  \code{\link[PBSsatellite]{print.ncdfData}},
  \code{\link[PBSsatellite]{read.ncdfData}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% use one of  RShowDoc("KEYWORDS")
\keyword{file}% __ONLY ONE__ keyword per line

