\name{km2lat}
\alias{km2lat}
\title{Convert long-lat to km-km}
\description{The function computes the latitude from given
  the distance from a reference point. See also \code{\link{km2lon}} and \code{\link{COn0E65N}}.}
\usage{km2lat(x, y, x.centre=0, y.centre=65)}
\arguments{
    \item{x}{not used.}
    \item{y}{distance from reference latitude in meridional direction.}
    \item{x.centre}{reference longitude.}
    \item{y.centre}{reference latitude.}
 }
\value{real}
\author{R.E. Benestad}}
\examples{
library(clim.pact)
data(oslo.t2m)
print(c(oslo.t2m$lon,oslo.t2m$lat))
#[1] 10.71667 59.95000
xy<-COn0E65N(oslo.t2m$lon,oslo.t2m$lat)
oslo.t2m$lon<-xy$x
oslo.t2m$lat<-xy$y
print(c(oslo.t2m$lon,oslo.t2m$lat))
#[1]  595.4086 -560.3004
lon<-km2lon(oslo.t2m$lon,oslo.t2m$lat,x.centre=0,y.centre=65)
lat<-km2lat(oslo.t2m$lon,oslo.t2m$lat,x.centre=0,y.centre=65)
print(c(lon,lat))
#[1] 10.71667 59.95000
}
\keyword{manip}
