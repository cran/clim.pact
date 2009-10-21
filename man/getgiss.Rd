\name{getgiss}
\alias{getgiss}
\alias{getnarp}
\title{Retrieve station record from the GISS or NARP data set from URL.}
\description{Reads the data from GISS available at URL:
  \url{http://www.giss.nasa.gov/data/update/gistemp/} or NARP data 
  from \url{http://projects.met.no/~narp/data_index.html} (\url{http://thule.oulu.fi/narp/}).  
Also see \code{\link{getnordklim}}, \code{\link{getdnmi}} and \code{\link{getnacd}}.

Use eg \code{options(timeout=120)} if problems with retreaving the
GISS data. Some times, the code doesn't manage to read the data, and
returns with the message: "cannot open: HTTP status was '404 Not
Found'". Try again, since slow connections cause this problem.
One trick that often helps is to use the function to get the name of
the station (\code{stations <- getgiss()}), then use a browser
(Mozilla) to open the data file from
\url{http://www.giss.nasa.gov/data/update/gistemp/station_data/}, and
then call \code{getgiss} again with the station details. It is easier to access
the NARP data.
}
\usage{
getgiss(stnr=NULL,location=NULL,lon=NULL,lat=NULL,stations=NULL,
        silent=FALSE)
getnarp(stnr=NULL,location=NULL,lon=NULL,lat=NULL,stations=NULL,
        silent=TRUE,ele=101)
}

\arguments{
  \item{stnr}{Station number.}
  \item{location}{Name of locvation.}
  \item{lon}{longitude. If both \code{lon} and \code{lat} are
  specified and the other arguments are not, then find the closest station.}
  \item{lat}{latitude.}
  \item{stations}{a list object with details of GISS stations
  (Prescribing this saves time).}
  \item{silent}{For verbose use.}
  \item{ele}{Element code: see \url{http://projects.met.no/~narp/data_index.html}.}
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
\author{R.E. Benestad}
\examples{
\dontrun{
#GISS:
obs.oxford <- getgiss(location="Oxford")  # Takes longer time
stations <- getgiss()
obs.oxford <- getgiss(location="Oxford",stations=stations) #Quicker
plotStation(obs.oxford)
obs.broome <- getgiss(stnr="501942030004",stations=stations)
obs.120E.40S <- getgiss(lon=120,lat=-40,stations=stations)

#NARP:
> obs <- getnarp(4360)
#[1] "Retrieving the data from URL http://projects.met.no/~narp/narp"
#[1] "Please be patient"
#[1] "Found Tasilaq stnr= 4360  lon= -37.63  lat= 65.6  country= G"
> obs <- getnarp(lon=-37,lat=65)
#[1] "Retrieving the data from URL http://projects.met.no/~narp/narp"
#[1] "Please be patient"
#[1] "Find the nearest station to -37E and 65N."
#[1] "Found Tasilaq stnr= 4360  lon= -37.63  lat= 65.6  country= G"
> obs <- getnarp(lon=10,lat=60)
#[1] "Retrieving the data from URL http://projects.met.no/~narp/narp"
#[1] "Please be patient"
#[1] "Find the nearest station to 10E and 60N."
#[1] "Found Oslo stnr= 18700  lon= 10.72  lat= 59.95  country= N"
}
}
\keyword{file}
