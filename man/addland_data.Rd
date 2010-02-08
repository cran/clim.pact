\name{addland1}
\alias{continental coast line}
\alias{lat.cont}
\alias{lon.cont}
\alias{addland2}
\title{Coordinates of coast line.}
\description{Coordinates of coast line: \code{lat.cont} and
  \code{lon.cont} give the latitude and longitude coordinates of the
  coasts (continents).}
\usage{
data(addland1)
data(addland2)}
\format{Two vectors, lat.cont and lon.cont, of length 124622.}
\source{Ferret web site}
\references{URL \url{http://ferret.wrc.noaa.gov/Ferret/}}
\examples{
library(clim.pact)
data(addland1)
ls()                 # [1] "lat.cont" "lon.cont"
}
\keyword{datasets}
