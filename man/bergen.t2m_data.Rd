\name{bergen.t2m}
\alias{bergen.t2m}
\title{Monthly mean temperature in Bergen.}
\description{A station record of monthly mean temperature from Bergen-Florida.}
\usage{data(bergen.t2m)}
\format{a <- list of "monthly.station.record" class:
    \tabular{ll}{
    val \tab The monthly values (a 12-column matrix with one
    column for each year). \cr
    station \tab station number. \cr
    yy \tab The years of observation (vector). \cr
    lat,lon \tab Latitude and longitude of the location. \cr
    x.0E65N,y.0E65N \tab Distance in km from 0E, 65N. \cr
    location \tab Name of location . \cr
    wmo.no \tab WMO number. \cr
    start \tab Start of observatins from this location. \cr
    yy0 \tab First year in current record. \cr
    ele \tab Code of theelement. \cr
    obs.name \tab Name of the element. \cr
    unit \tab Unit of the element. \cr
    country \tab The country in which the location is
    located. \cr
    quality \tab Code/description for data quality. \cr
    found \tab Flag: T - the data requested was found. \cr
    ref \tab Reference for the data set. \cr
  }}
\source{The Nordklim data set: \url{http://www.smhi.se/hfa_coord/nordklim/}}
\references{Tuomenvirta et al. (2001), "Nordklim data set 1.0", DNMI
  KLIMA 08/01, pp. 26; \cr
  The Norwegian Meteorological Institute, P.O. Box 43, 0313
Oslo, Norway (http://www.met.no).}
\examples{
}
\keyword{datasets}