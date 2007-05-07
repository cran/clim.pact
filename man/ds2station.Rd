\name{ds2station}
\alias{ds2station}
\title{Convert ds to station object}
\description{Convert 'ds' object to 'station' object so that
downscaled scenarios can be handled like station series.} 
\usage{ds2station(x,what="scenario")}
\arguments{
    \item{x}{a 'ds' object (see \code{\link{DS}}).}
    \item{what}{either 'scenario' or 'calibration'}
  }
\value{a station object}
\author{R.E. Benestad}
\examples{
data(eof.c)
data(oslo.t2m)
ds <- DS(preds=eof.c, oslo.t2m,plot=FALSE)
oslo.dm <- ds2station(ds)
plotStation(oslo.dm)
}
\keyword{manip}
