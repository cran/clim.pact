\name{plotField}
\alias{plotField}
\title{plotField}
\description{Produce 2D plots like maps and Hovmuller diagrams. A
  poor man's version of Ferret's
  \url{http://ferret.wrc.noaa.gov/Ferret/} plot
  function. \code{plot.field} is a high level command that utilises
  \code{\link{mapField}} or \code{\link{grd.box.ts}} whenever appropriate.}  
\usage{plotField(x,lon=NULL,lat=NULL,tim=NULL,mon=NULL,
                       col="black",lty=1,lwd=1,what="ano")}
\arguments{
  \item{x}{A field object.}
  \item{lon}{Position for longitude. One of \code{lon, lat, tim} must be set.}
  \item{lat}{Position for latitude.}
  \item{tim}{Position for time.}
  \item{mon}{Month to extract.}
  \item{col}{Contour line colour.}
  \item{lty}{Contour line type.}
  \item{lwd}{Contour line width.}
  \item{what}{Choose between "ano" - anomaly; "abs" absolute; "cli"
    climatological.} 
}
\value{}
\author{R.E. Benestad}}
\examples{
}
\keyword{hplot}