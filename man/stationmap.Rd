\name{stationmap}
\alias{stationmap}
\title{Plot climate station map.}
\description{Plot a map showing the locations of NACD (\code{\link{getnacd}})
  and Nordklim (\code{\link{getnordklim}}) station network (monthly data).}
\usage{stationmap(ele=101,NORDKLIM=TRUE,NACD=TRUE,silent=TRUE)}
\arguments{
  \item{ele}{The code for element (101 = T(2m), 401=SLP, 601= precip).}
  \item{NORDKLIM}{Show Nordklim stations.}
  \item{NACD}{Show NACD stations.}
  \item{silent}{Print station names if false.}
\value{}
\author{R.E. Benestad}}
\examples{
\dontrun{
stationmap()
}
}
\keyword{hplot}
