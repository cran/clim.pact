\name{plotField}
\alias{plotField}
\title{plotField}
\description{Produce 2D plots like maps and Hovmuller diagrams. A
  poor man's version of Ferret's
  \url{http://ferret.wrc.noaa.gov/Ferret/} plot
  function. \code{plot.field} is a high level command that utilises
  \code{\link{mapField}} or \code{\link{grd.box.ts}} whenever appropriate.}  
\usage{plotField(x,lon=NULL,lat=NULL,tim=NULL,mon=NULL,val.rng=NULL,
                       col="black",col.coast="grey",lty=1,lwd=1,what="ano",
                       type="s",pch=26,my.col=NULL,add=FALSE)}
\arguments{
  \item{x}{A field object.}
  \item{lon}{Position for longitude. One of \code{lon, lat, tim} must be set.}
  \item{lat}{Position for latitude.}
  \item{tim}{Position for time as time index (1:length(tim).}
  \item{mon}{Month to extract.}
  \item{val.rng}{Valid range: interval used for colour scale.}
  \item{col}{Contour line colour for levels.}
  \item{col.coast}{Contour line colour for coast lines.}
  \item{lty}{Contour line type.}
  \item{lwd}{Contour line width.}
  \item{what}{Choose between "ano" - anomaly; "abs" absolute; "cli"
    climatological.}
  \item{type}{same as in plot(), only works for single time series.}
  \item{pch}{same as in plot(), only works for single time series.}
  \item{my.col}{colour palette (see \code{link{rgb}}).}
  \item{add}{TRUE adds a time series to previous  (see \code{link{grd.box.ts}}).}
}
\value{}
\author{R.E. Benestad}}
\examples{
\dontrun{
skt <- retrieve.nc("skt.mon.mean.nc",x.rng=c(-90,50),y.rng=c(0,75))

# Maps of monthly mean skin temperatures:
plotField(skt,tim=1,val.rng=c(-20,20))
dev2bitmap("ncep.skt_194801.jpg",type="jpeg")


plotField(skt,tim=100,col="blue",col.coast="darkgreen",val.rng=c(-10,10))

# For adding extra points/contours:

    # From filled.contour in base
    mar.orig <- (par.orig <- par(c("mar","las","mfrow")))$mar
    on.exit(par(par.orig))

    w <- (3 + mar.orig[2]) * par('csi') * 2.54
    layout(matrix(c(2, 1), nc=2), widths=c(1, lcm(w)))
   
    par(las = 1)
    mar <- mar.orig
    mar[4] <- 1
    par(mar=mar)
# End of section affecting the window set up.

points(0,50,pch=21,col="red")
grid()
dev2bitmap("ncep.skt_195604.jpg",type="jpeg")


# A hovmuller diagram:
plotField(skt,lon=0,val.rng=c(-10,10))
dev2bitmap("ncep.skt_lontim.jpg",type="jpeg")


# A single time series:
plotField(skt,lon=-20,lat=50)

}
}
\keyword{hplot}
