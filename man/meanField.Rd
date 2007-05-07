\name{meanField}
\alias{meanField}
\title{Mean field}
\description{Computes a map showing the mean field values.}
\usage{meanField(x,lon.rng=NULL,lat.rng=NULL,t.rng=NULL,mon=NULL)}
\arguments{
    \item{x}{A field object.}
    \item{lon.rng}{Extract the longitude interval.}
    \item{lat.rng}{Extract the latitude interval.}
    \item{t.rng}{Extract the time interval. Can be numerical (year) or character (date, e.g. "01-Jan-1998"}
    \item{mon}{Month (1-12) or season (1-12) to extract if !is.null.}
 }
\value{A map object}
\author{R.E. Benestad}
\examples{
\dontrun{
slp <- retrieve.nc("ncep_slp.nc",x.rng=c(5,12),y.rng=c(58,63))
mslp <- meanField(slp)
}
}
\keyword{manip}
