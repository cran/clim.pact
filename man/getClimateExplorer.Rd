\name{getClimateExplorer}
\alias{getClimateExplorer}
\title{Retrieve station record through the KNMI ClimateExplorer.}
\description{Retrieve data through: URL \url{http://climexp.knmi.nl/}. Also see
  \code{\link{getdnmi}}, \code{\link{getnordklim}}, and \code{\link{getnacd}}.
}
\usage{
getClimateExplorer(location,URL="http://climexp.knmi.nl/data/",ele.c="t")
}

\arguments{
  \item{location}{ClimateExplorer station number.}
  \item{ele.c}{name of element [e.g. avail.elem(), or 't', 'p', 'l'].}
  \item{URL}{URL for the data report.}
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
  }}
\author{R.E. Benestad}
\examples{
\dontrun{
greenwich <- getClimateExplorer(3779.1)
plotStation(greenwich,what="t")
all.station.list <- getClimateExplorer()
}
}
\keyword{file}
