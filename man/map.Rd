\name{map}
\alias{map}
\title{Produce a map}
\description{Produces maps.}
\usage{map(x,y=NULL,col="black",lwd=1,lty=1,sym=TRUE,
               plot=TRUE,inv.col=FALSE)}
\arguments{
    \item{x}{A map object.}
    \item{y}{A map object. If given, \code{map.map} plots the
      difference: x - y}
    \item{col}{Colour of contours.}
    \item{lwd}{Contour line width.}
    \item{lty}{Contour line style.}
    \item{sym}{Symmetry: if True, use zlimits c(-|max|,+|max|).}
    \item{plot}{TRUE gives graphics}
    \item{inv.col}{Inverse color scheme (e.g. 'TRUE' gives red for drier
    and blue for wetter conditions)}.
 }
\value{A map object}
\author{R.E. Benestad}}
\examples{
}
\keyword{hplot}