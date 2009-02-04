\name{grd.box.ts}
\alias{grd.box.ts}
\title{Grid box time series}
\description{Plots time series from field objects, eg from
  \code{\link{retrieve.nc}}. Uses a bilinear interpolation scheme
  (\code{link{interp}}) from the \code{akima} package.}
\usage{
grd.box.ts(x,lon,lat,what="abs",greenwich=TRUE,mon=NULL,
                       col="grey10",lwd=1,lty=1,pch=".",add=FALSE,
                       filter=NULL,type="l",main=NULL,sub=NULL,xlab=NULL,
                       ylab=NULL,xlim=NULL,ylim=NULL)
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
    \item{type}{same as in plot(), only works for single time series.}
    \item{main}{Preset main title (see \code{\link{plot}}).}
    \item{sub}{subtitle}
    \item{xlab}{x label}
    \item{ylab}{y label}
    \item{xlim}{see \code{\link{plot}}}
    \item{ylim}{see \code{\link{plot}}}

}
\value{A station object with interpolated values (see \code{\link{station.obj}} and
  \code{\link{station.obj.dm}}) that can be plotted using eg \code{\link{plotStation}}.    

  }
\author{R.E. Benestad}
\examples{
\dontrun{
slp <- retrieve.nc("ncep_slp.nc")
grd.box.ts(slp,0,60,what="ano",mon=1)
}
}
\keyword{ts}





