\name{anomaly.field}
\alias{anomaly.field}
\title{Anomalies of a field object.}
\description{Estimates the anomalies of a field object.  Also see \code{\link{anomaly.field}}.}
\usage{anomaly.field(x,period=NULL)}
\arguments{
  \item{x}{A field object.}
  \item{period}{Period to use as climatology 
           NULL -> the entire series as reference clim.}
 }
\value{A field object.}
\author{R.E. Benestad}}
\examples{
\dontrun{
slp<-retrieve.nc("ncep_slp.nc")
slp.a<-anomaly.field(slp)
}
}
\keyword{hplot}