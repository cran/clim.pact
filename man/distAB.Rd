\name{distAB}
\alias{distAB}
\title{Distance between two points on Earth}
\description{The function returns the distance between two points on
  Earth given by the lon-lat coordinates. The distance is computed using
the formula: d = a * theta, and |x1|*|x2| * cos(theta) = inner-product(x1,x2).}
\usage{distAB(lon,lat,lons,lats,a=6.378e06)}
\arguments{
    \item{lon}{Longitude of reference point (degrees East).}
    \item{lat}{Latitude of reference point (degrees North).}
    \item{lons}{Longitude of points of interest (vector) (degrees East).}
    \item{lats}{Latitude of points of interest (vector) (degrees
      North).}
    \item{a}{Radius of the Earth.}
 }
\value{A real value: units=meters.}
\author{R.E. Benestad}}
\examples{
distAB(10,60,5,58)    # [1] 362802.3
distAB(0,0,180,0)     # [1] 20037078
distAB(0,90,0,-90)    # [1] 20037078
}
\keyword{math}