\name{map.eof}
\alias{map.eof}
\title{Map eof}
\description{Draws maps of the spatial structures described by the EOF
  (\code{\link{EOF}} patterns. Is similar to \code{\link{plotEOF}}, but
  only plots the spatial information. Useful for comparing the spatial
  patterns in different EOFs.} 
\usage{map.eof(x,i.eof=1,nlevs=5,add=FALSE,
            col=c("red","blue","darkgreen","steelblue"),lwd=2,lty=1)}
\arguments{
    \item{x}{ A field object. }
    \item{i.eof}{The EOF to plot.}
    \item{nlevs}{Number of contour levels. }
    \item{add}{Add a map on pre-existing map - see contour.}
    \item{col}{Colour.}
    \item{lwd}{Line width.}
    \item{lty}{Line type.}
  }
\value{A map object (currently, only the last field in a mixed.field object).}
\author{R.E. Benestad}}
\examples{
library(clim.pact)
data(eof.slp)
map.eof(eof.slp)
map.eof(eof.slp,i.eof=2,col="blue",add=TRUE)
}
\keyword{hplot}
