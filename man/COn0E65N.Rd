\name{COn0E65N}
\alias{COn0E65N}
\alias{lat}
\alias{lon}
\title{Convert long-lat to km-km}
\description{The function returns the distance in km from 0E 65N given
  the longitude and latitude. See also \code{\link{km2lon}} and \code{\link{km2lat}}.}
\usage{COn0E65N(lon, lat)}
\arguments{
 \item{lon}{longitude}
 \item{lat}{latitude}
}
\value{list(y=latitudes distance,x= longitudes distance)}
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
