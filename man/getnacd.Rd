\name{getnacd}
\alias{getnacd}
\alias{nacd.meta}
\alias{meta.nacd}
\title{Retreave station record from the NACD set.}
\description{Retrieve station record from the North Atlantic Climate
  Data (NACD) set: URL \url{http://www.dmi.dk/}. Also see
  \code{\link{getdnmi}} and \code{\link{getnordklim}}.
}
\usage{
getnacd(location="prompt",ele.c='101',ascii=FALSE,silent=FALSE,direc="data")
}

\arguments{
  \item{location}{name of climate station location.}
  \item{ele.c}{name of element [e.g. avail.elem(), or 't2m', 'rr', 'slp'].}
  \item{ascii}{Flag. T -> force ascii read, otherwise look for
    R-formatted verion (faster).}
  \item{silent}{Flag. F -> print error messages.}
  \item{direc}{name of directory in which the data are stored.}
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
\author{R.E. Benestad}}
\examples{
\dontrun{
helsinki.rr <- getnacd("helsinki",ele=601)
obs.t2m <- getnacd()
}
}
\keyword{file}