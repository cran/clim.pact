\name{grd.box.ts}
\alias{grd.box.ts}
\title{Grid box time series}
\description{Plots time series from field objects, eg from \code{\link{retrieve.nc}}.}
\usage{
grd.box.ts(x,lon,lat,what="abs",greenwich=TRUE,mon=NULL,
                       col="grey10",lwd=1,lty=1,pch=26,add=FALSE,
                       filter=NULL)
}
\arguments{
    \item{x}{A field object.}
    \item{lon}{Longitude to plot.}
    \item{lat}{Latitude to plot.}
    \item{what}{What to draw: "ano"-> anomalies, "cli"-> climatological
      values, "abs" -> absolute values.}
    \item{greenwich}{Maps centre on the Greenwich meridian.}
    \item{mon}{Month to extract}
    \item{col}{Colour.}
    \item{lwd}{Line width}
    \item{lty}{Line style.}
    \item{pch}{Plot character.}
    \item{add}{'TRUE' adds curve to old plot.}
    \item{filter}{If not NULL, this is a vector specifying window
      weighting argument of same name in \code{link{filter}}.}
}
\value{
    
    \tabular{ll}{
      avail.preds \tab vector of characters \cr
      avail.locs \tab a list with name, lons, lats,
      country, ident \cr
      avail.elem \tab a list with data.set, ele, name \cr}}
\author{R.E. Benestad}}
\examples{
\dontrun{
slp <- retrieve.nc("ncep_slp.nc")
grd.box.ts(slp,0,60,what="ano",mon=1)
}
}
\keyword{ts}
