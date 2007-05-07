\name{stationmap}
\alias{stationmap}
\title{Plot climate station map.}
\description{Plot a map showing the locations of NACD (\code{\link{getnacd}})
  and Nordklim (\code{\link{getnordklim}}) station network (monthly data).}
\usage{stationmap(ele=101,NORDKLIM=TRUE,NACD=TRUE,silent=TRUE,names=FALSE,
                  name.len=4,x.offs=0.1,y.offs=-0.5,str.cex=0.7,
                  countries=NULL,x.rng=NULL,y.rng=NULL)}
\arguments{
  \item{ele}{The code for element (101 = T(2m), 401=SLP, 601= precip).}
  \item{NORDKLIM}{Show Nordklim stations.}
  \item{NACD}{Show NACD stations.}
  \item{silent}{Print station names if false.}
  \item{names}{Plot station names if true.}
  \item{name.len}{Station name text length.}
  \item{x.offs}{Station name text position.}
  \item{y.offs}{Station name text position.}
  \item{str.cex}{Station name text size.}
  \item{countries}{Countries to include.}
  \item{x.rng}{Longitudes to include.}
  \item{y.rng}{Latitudes to include.}
  }
\value{}
\author{R.E. Benestad}
\examples{
\dontrun{
stationmap()
}
}
\keyword{hplot}
