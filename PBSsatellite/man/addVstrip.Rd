\name{addVstrip}
\alias{addVstrip}
\title{
  Add Vertical Strip
}
\description{
  Add a vertical colour strip as a legend.
}
\usage{
addVstrip(x, y, zlim, col, title, xwidth=0.02, yheight=0.5)
}
\arguments{
  \item{x}{\code{numeric} -- relative x-coordinate to place the strip.}
  \item{y}{\code{numeric} -- relative y-coordinate to place the strip.}
  \item{zlim}{\code{numeric} -- limits of the z-value used to label the colour strip.}
  \item{col}{\code{character} -- vector of colours to create the colour strip.}
  \item{title}{\code{character} -- title to appear abouve the strip.}
  \item{xwidth}{\code{numeric} -- relative width of strip.}
  \item{yheight}{\code{numeric} -- relative height of strip.}
}
\details{
  Primitive function; for improved version see \code{PBStools::addStrip}.
}
\value{
  List of strip-related measurements.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2020-11-05}
}
\seealso{
In package \pkg{graphics}:
  \code{\link[graphics]{polygon}}\cr
In package \pkg{PBStools}:
  \code{\link[PBStools]{addStrip}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{aplot}
\keyword{color}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{legend}
% \concept{ ~cpt2 }
% Use only one concept per line.
