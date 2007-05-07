\name{mapField}
\alias{mapField}
\title{MapField}
\description{Draws maps of fields in a field object, eg read using \code{\link{retrieve.nc}}.}
\usage{mapField(x,l=NULL,greenwich=TRUE,plot=TRUE,
                     what="ano",method="nice",val.rng=NULL,
                     col="black",col.coast="grey",lwd=2,lty=1,
                     add=FALSE,las = 1,levels=NULL,xlim=NULL,ylim=NULL,newFig=TRUE)
}
\arguments{
    \item{x}{ A field object. }
    \item{l}{The field to map. Default: the last field in the record.}
    \item{what}{What to draw: "ano"-> anomalies, "cli"-> climatological
    values, "abs" -> absolute values. }
  \item{method}{"nice" -> \code{filled.contour}, otherwise use
    \code{image}.}
    \item{val.rng}{Valid range: interval used for colour scale.}
    \item{greenwich}{'TRUE' centres on Greenwich meridian.}
    \item{plot}{'TRUE' produces graphic.}
    \item{col}{Contour line colour for levels.}
    \item{col.coast}{Contour line colour for coast lines.}
    \item{lwd}{Line width.}
    \item{lty}{Line type.}
    \item{add}{Adds map to old figure.}
    \item{las}{See \code{\link{par}}}
    \item{levels}{See \code{\link{contour}}}
    \item{xlim}{see \code{\link{plot}}}
    \item{ylim}{see \code{\link{plot}}}
    \item{newFig}{TRUE call \code{\link{newFig}} for new figure.}
  }
\value{}
\author{R.E. Benestad}
\examples{
library(clim.pact)
\dontrun{
skt<-retrieve.nc("skt.mon.mean.nc",
                 x.rng=c(-90,50),y.rng=c(0,75))
bitmap("ncep.skt.jpg",type="jpeg")
mapField(skt)
dev.off()
}
}
\keyword{hplot}
