\name{plotStation}
\alias{plotStation}
\title{Plots monthly station records.}
\description{Plots data in the monthly station records. The data
  may be read through \link{getnacd}, \link{getnordklim},
  \link{getdnmi}, or created using \link{station.obj}. The commands
  \link{avail.elem}, \link{avail.locs} can be used to identify the
  station records available (in a given subdirectory).
}
\usage{plotStation(obs,l.anom=TRUE,mon=NULL,
                        leps=FALSE,out.dir="output",what="b",
                        type="l",pch=26,col="black",lwd=3,lty=3,add=FALSE)}

\arguments{
  \item{obs}{A climate station series. }
  \item{l.anom}{flag: T -> plot anomalies. }
  \item{mon}{select month to plot, A value of 0 plots all
    months. mon=c(12,1,2) plots the DJF mean.  }
  \item{leps}{Flag: T -> produce EPS files (hard copy). }
  \item{out.dir}{Directory where to store hard copies.}
  \item{what}{"t"=timeseries; "d"=distribution; "b"=both.}
  \item{type}{same as in \code{link{plot}}.}
  \item{pch}{same as in \code{link{plot}}.}
  \item{col}{colour, same as in \code{link{plot}}.}
  \item{lwd}{line width, same as in \code{link{plot}}.}
  \item{lty}{line type, same as in \code{link{plot}}.}
  \item{add}{TRUE adds a time series to previous  (see \code{link{grd.box.ts}}).}
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
data(tromsoe.t2m)
plotStation(tromsoe.t2m)
}
}
\keyword{hplot}
