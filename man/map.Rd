\name{map}
\alias{map}
\title{Produce a map}
\description{Produces maps.}
\usage{map(x,y=NULL,col="black",lwd=1,lty=1,sym=TRUE,
               plot=TRUE,inv.col=FALSE,add=FALSE,las = 1,
                    levels=NULL,main=NULL,sub=NULL,xlim=NULL,ylim=NULL)}
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
    \item{add}{Add to plot: add contours only.}
    \item{las}{see \code{\link{par}}}
    \item{levels}{see \code{\link{contour}}}
    \item{main}{see \code{\link{plot}}}
    \item{sub}{see \code{\link{plot}}}
    \item{xlim}{see \code{\link{plot}}}
    \item{ylim}{see \code{\link{plot}}}
 }
\value{A map object}
\author{R.E. Benestad}}
\examples{
}
\keyword{hplot}
