\name{station.obj.dm}
\alias{station.obj.dm}
\alias{daily.station.record}
\title{Make daily climate station series object.}
\description{Create a station object for use as predictand in empirical
  downscaling on monthly data. Also see \code{\link{station.obj}}.}
\usage{station.obj.dm(t2m,precip,dd,mm,yy,
                       obs.name=NULL,unit=NULL,ele=NULL,
                       station=NULL,lat=NULL,lon=NULL,alt=NULL,
                       location="unspecified",wmo.no=NULL,
                       start=NULL,yy0=NULL,country=NULL,ref=NULL)}
\arguments{
  \item{t2m}{a vector holding daily mean temperature. }
  \item{precip}{a vector holding daily precipitation. }
  \item{dd}{a vector holding day of month.}
  \item{mm}{a vector holding the month.}
  \item{yy}{a vector holding the year.}
  \item{obs.name}{the name of observation: 
  eg c("Daily mean temperature","Daily precipitation"). }
  \item{unit}{the unite of observation: eg c("deg C","mm/day").}
  \item{ele}{element code: eg c("tam","rr").}
  \item{station}{local (national) station number.}
  \item{lat}{latitude.}
  \item{lon}{longitude.}
  \item{alt}{altitude.}
  \item{location}{name of location. }
  \item{wmo.no}{WMO number of station.}
  \item{start}{start of measurements.}
  \item{yy0}{first year of record.}
  \item{country}{name of country.}
  \item{ref}{reference to the data.}
  }
  \value{a "daily.station.record"-class object.
  \tabular{ll}{
  t2m \tab a vector holding daily mean temperature. \cr
  precip \tab a vector holding daily precipitation. \cr
  day \tab a vector holding day of month.\cr
  month \tab a vector holding the month.\cr
  year \tab a vector holding the year.\cr
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
  ref \tab reference to the data.\cr}
  }
\author{R.E. Benestad}}
\examples{
\dontrun{
blindern.raw <-read.table("~/data/stations/blindern_rr_day.dat",header=TRUE)
blindern.raw$rr[blindern.raw$rr < 0] <- NA
yy <- floor(blindern.raw$yyyymmdd/10000)
mm <- floor(blindern.raw$yyyymmdd/100) - 10000*yy
dd <- blindern.raw$yyyymmdd - 100*mm  - 10000*yy
blindern <- station.obj.dm(t2m=rep(NA,length(blindern.raw$rr)),
                           precip=blindern.raw$rr,
                           dd=dd,mm=mm,yy=yy,
                           obs.name=c("T(2m)","recip"),
                           unit=c("deg C","mm/day"),ele=NULL,
                           station=18700,lat=59.95,lon=10.71,alt=94,
                           location="Oslo-Blindern",wmo.no=NULL,
                           start=NULL,yy0=1937,country="Norway",
                           ref="www.met.no")
}
}
\keyword{data}
