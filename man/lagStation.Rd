\name{lagStation}
\alias{lagStation}
\title{Introduce a lag in station values}
\description{Shifts the station values in time: introduce a lag} 
\usage{lagStation(x,lag=0)}
\arguments{
    \item{x}{A station object (see \code{\link{station.obj}}).}
    \item{lag}{lag}
  }
\value{Astation object}
\author{R.E. Benestad}}
\examples{
data(oslo.t2m)
oslo.t2m.1 <- lagStation(oslo.t2m,1)
plotStation(oslo.t2m,mon=1,what="t")
plotStation(oslo.t2m.1,add=TRUE,mon=1,what="t",col="darkblue")
}
\keyword{manip}
