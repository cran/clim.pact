\name{anomaly.field}
\alias{anomaly.field}
\alias{daily2monthly.field}
\title{Anomalies of a field object.}
\description{Estimates the anomalies of a field object.  Also see \code{\link{anomaly.field}}.}
\usage{
anomaly.field(x,period=NULL)
daily2monthly.field(field,min.days.month=20,method="colMeans",na.rm=TRUE)
}
\arguments{
  \item{x}{A field object.}
  \item{period}{Period to use as climatology 
    NULL -> the entire series as reference clim.}
  \item{field}{}
  \item{min.days.month}{}
  \item{method}{}
  \item{na.rm}{}
 }
\value{A field object.}
\author{R.E. Benestad}
\examples{
\dontrun{
slp<-retrieve.nc("ncep_slp.nc")
slp.a<-anomaly.field(slp)
}
}
\keyword{hplot}
