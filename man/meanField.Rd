\name{meanField}
\alias{optimal interpolation}
\alias{meanField}
\title{Optimal interpolation}
\description{.}
\usage{meanField(x,lon.rng=NULL,lat.rng=NULL,t.rng=NULL)}
\arguments{
    \item{x}{A field object.}
    \item{lon.rng}{Extract the longitude interval.}
    \item{lat.rng}{Extract the latitude interval.}
    \item{t.rng}{Extract the time interval.}
 }
\value{A map object}
\author{R.E. Benestad}}
\examples{
\dontrun{
slp <- retrieve.nc("ncep_slp.nc",x.rng=c(5,12),y.rng=c(58,63))
mslp <- meanField(slp)
}
}
\keyword{manip}