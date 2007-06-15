\name{getnordklim}
\alias{getnordklim}
\alias{meta}
\alias{getHadObs}
\title{Retrieve station record from the Nordklima set or HadObs.}
\description{getnordklim reads the data from Nordklim available at URL:
  \url{http://www.smhi.se/hfa_coord/nordklim/}.  Also see
  \code{\link{getdnmi}} and \code{\link{getnacd}}.

  
}
\usage{
getnordklim(location=NULL,ele.c='101',ascii=FALSE,silent=FALSE,direc="data")
getHadObs(what="CET")
}

\arguments{
  \item{location}{ name of climate station location. NULL returns the
    list of all names.}
  \item{ele.c}{ name of element [e.g. avail.elem(), or 't2m', 'rr', 'slp'].}
  \item{ascii}{ Flag. T -> force ascii read, otherwise look for
    R-formatted verion (faster).}
  \item{silent}{ Flag. F -> print error messages.}
  \item{direc}{ name of directory in which the data are stored.}
  \item{what}{"CET": reads the central England Temperature; "EWP": reads Monthly England and Wales precipitation}
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
helsinki.rr <- getnordklim("helsinki",ele=601)
}
}
\keyword{file}
