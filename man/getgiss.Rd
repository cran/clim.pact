\name{getgiss}
\alias{getgiss}
\title{Retrieve station record from the GISS data set from URL.}
\description{Reads the data from GISS available at URL:
  \url{http://www.smhi.se/hfa_coord/nordklim/}.  Also see
  \code{\link{getdnmi}} and \code{\link{getnacd}}.

Use eg \code{options(timeout=120)} if problems with retreaving the
data. Some times, the code doesn't manage to read the data, and
returns with the message: "cannot open: HTTP status was '404 Not
Found'". Try again, since slow connections cause this problem.
One trick that often helps is to use the function to get the name of
the station (\code{stations <- getgiss()}), then use a browser
(Mozilla) to open the data file from
\url{http://www.giss.nasa.gov/data/update/gistemp/station_data/}, and
then call \code{getgiss} again with the station details.
}
\usage{getgiss(stnr=NULL,location=NULL,lon=NULL,lat=NULL,stations=NULL,silent=FALSE)}

\arguments{
  \item{stnr}{Station number.}
  \item{location}{Name of locvation.}
  \item{lon}{longitude. If both \code{lon} and \code{lat} are
  specified and the other arguments are not, then find the closest station.}
  \item{lat}{latitude.}
  \item{stations}{a list object with details of GISS stations
  (Prescribing this saves time).}
  \item{silent}{For verbose use.}
}
  \value{a <- list of "monthly.station.record" class:
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
  }
or a list of station details if the requested station was not found.
}
\author{R.E. Benestad}}
\examples{
\dontrun{
obs.oxford <- getgiss(location="Oxford")  # Takes longer time
stations <- getgiss()
obs.oxford <- getgiss(location="Oxford",stations=stations) #Quicker
plotStation(obs.oxford)
obs.broome <- getgiss(stnr="501942030004",stations=stations)
obs.120E.40S <- getgiss(lon=120,lat=-40,stations=stations)
}
}
\keyword{file}
