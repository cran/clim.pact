\name{anomaly.station}
\alias{anomaly.station}
\title{Anomaly.station}
\description{Computes anomalies of a station series by subtracting the
  climatology. The climatology is estimated either by taking the average
  of the respective months over a given reference period or a least
  squares fit to the 6 leading harmonics, depending on the
  appropriateness. Also see \code{\link{anomaly.field}}.}
\usage{
anomaly.station(obs,period=c(1961,1990),param=c("t2m","precip"))
}
\arguments{
    \item{obs}{Monthly station series}
    \item{period}{Period to use as climatology 
           NULL -> the entire series as reference clim.}
    \item{param}{For daily objects: identifies the names of the parameters.}
  }
\value{Station series object.}
\author{R.E. Benestad}}
\examples{
library(clim.pact)
data(oslo.t2m)
oslo.t2ma<-anomaly.station(oslo.t2m)
}
\keyword{manip}
