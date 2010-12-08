\name{bergen.dm}
\alias{bergen.dm}
\title{Daily Bergen record.}
\description{A station record of daily mean temperature and daily
  precipitation from Bergen-Florida.}
\usage{data(bergen.dm)}
\format{a "daily.station.record"-class object.
  \tabular{ll}{
  t2m \tab a vector holding daily mean temperature. \cr
  precip \tab a vector holding daily precipitation. \cr
  dd \tab a vector holding day of month.\cr
  mm \tab a vector holding the month.\cr
  yy \tab a vector holding the year.\cr
  obs.name \tab the name of observation: \cr
  \tab eg c("Daily mean temperature","Daily precipitation"). \cr
  unit \tab the unite of observation: eg c("deg C","mm/day").\cr
  ele \tab element code: eg c("tam","rr").\cr
  station \tab local (national) station number.\cr
  lat \tab latitude.\cr
  lon \tab longitude.\cr
  alt \tab altitude.\cr
  location \tab name of location. \cr
  wmo.no \tab WMO number of station.\cr
  start \tab start of measurements.\cr
  yy0 \tab first year of record.\cr
  country \tab name of country.\cr
  ref \tab reference to the data.\cr}}
\source{The Norwegian Meteorological Institute, Climatology deivision.}
\references{The Norwegian Meteorological Institute, P.O. Box 43, 0313
Oslo, Norway (http://www.met.no).}
\examples{
data(bergen.dm)
}
\keyword{datasets}
