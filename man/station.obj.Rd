\name{station.obj}
\alias{station.obj}
\alias{monthly.station.record}
\alias{dat}
\title{Make monthly climate station series object.}
\description{Create a station object for use as predictand in empirical
  downscaling on monthly data.  Also see \code{\link{station.obj.dm}}.}
\usage{station.obj(x,yy,obs.name,unit,ele=NULL, mm=NULL,
                   station=NULL,lat=NULL,lon=NULL,alt=NULL,
                   location="unspecified",wmo.no=NULL,
                   start=NULL,yy0=NULL,country=NULL,ref=NULL)}
\arguments{
  \item{x}{the data: a matrix of 12 columns holding the observations of each
  calendar month: column 1 holds January values, col 2 holds February, .. 
        col 12 holds December values.}
  \item{yy}{A vector holding the year of observation of the same length
    as each of the 12 columns. Or a vector with the same length as the data if 
    \code{mm} is given.}
  \item{mm}{a vector of months with the same length as the data (optional).}
  \item{obs.name}{the name of observation (e.g. "Temperature").}
  \item{unit}{the unite of observation (e.g. "deg C").}
  \item{ele}{element code (e.g. 101).}
  \item{station}{local (national) station number.}
  \item{lat}{latitude.}
  \item{lon}{longitude.}
  \item{alt}{altitude.}
  \item{location}{name of location.}
  \item{wmo.no}{WMO number of station.}
  \item{start}{start of measurements.}
  \item{yy0}{first year of record.}
  \item{country}{name of country.}
  \item{ref}{reference to the data.}
  }
  \value{a "monthly.station.record"-class object.
  \tabular{ll}{
    val \tab The monthly values (a 12-column matrix with one
    column for each year)\cr
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
    found \tab Flag: T-> the data requested was found. \cr
    ref \tab Reference for the data set. \cr
  }
\author{R.E. Benestad}}
\examples{
\dontrun{
a <- read.table("data/bjornholt.dat",
                  col.names=c("station","year","month","rr",
                    "tam","sam","sdm","uum","pom","tax","tan"))
obs <- station.obj(x=a$rr,yy=a$year,mm=a$month,
                   obs.name="Precipitation",unit="mm",ele=601,
                   lat=60.03,lon=10.41,alt=360,
                   station=a$station[1],location="Bjornholt",
                   country="Norway",ref="met.no Climate data base")
plot(obs,mon=11)
}
}
\keyword{data}