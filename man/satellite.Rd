\name{satellite}
\alias{polar stereographic}
\alias{stereogr}
\alias{satellite}
\alias{map2sphere}
\alias{satelliteOld}
\title{Satellite view / polar stereographic}
\description{Produces polar stereographic maps / satellite
  views. \code{stereogr} is simpler, newer, and faster then \code{satellite}.}
\usage{
stereogr(map.obj,NH=TRUE,lat.0=0,inv.col=FALSE,levels=NULL,sym=TRUE,dr=0.01,main=NULL)
satelliteOld(map.obj,col="black",lwd=2,lty=1,add=FALSE,
                      las = 1,lon.0=NULL,lat.0=NULL,method="normal",
                      ni=100,nj=100,n.nearest=4,max.dist=3,landdata=addland2)
map2sphere(z,X=seq(-180,180,by=1),Y=seq(-90,90,by=1),pal="rainbow",nlevs=21,breaks=NULL,quick=FALSE,add=FALSE)
satellite(map,lon.0=0,lat.0=0,pal="rainbow",nlevs=21,breaks=NULL,quick=FALSE,add=FALSE)
}
\arguments{
    \item{map.obj}{a map object (\code{\link{mapField}}).}
    \item{col}{contour colours.}
    \item{lwd}{contour line width.}
    \item{lty}{contour line style.}
    \item{add}{FLAG: true adds contour onto old plot.}
    \item{las}{see \code{\link{par}}.}
    \item{lon.0}{Reference longitude: centre of map; NULL selects automatically.}
    \item{lat.0}{Reference latitude: centre of map; NULL selects
      automatically.}
    \item{method}{A choice between "normal","polarstereo", and "distance".}
    \item{ni}{Number of grid points along x-axis in new grid.}
    \item{nj}{Number of grid points along y-axis in new grid.}
    \item{n.nearest}{Number of points to use in re-gridding.}
    \item{max.dist}{The maximum inter-point distance used for re-gridding.}
    \item{landdata}{"addland" gives hiher resolution for coast lines, but is slower.}
    \item{inv.col}{Inverse color scheme.}
    \item{NH}{Northern Hemisphere?.}
    \item{levels}{contour/shading levels.}
    \item{main}{Main title - see \code{\link{plot}}.}
    \item{sym}{symmetric colour scale?.}
    \item{dr}{used for masking interpolated values in data voids. Look
      of result may be sensitive to spatial grid resolution.}
    \item{z}{The gridded data.}
    \item{X}{x-coordinate of z}
    \item{Y}{y-coordinate of z}
    \item{map}{mpa object to be mapped.}
    \item{pal}{palette.}
    \item{quick}{Faster high-resolution plot.}
    \item{nlevs}{number of levels for contour/shading levels}
    \item{breaks}{break points for contour/shading levels}
 }
\value{}
\author{R.E. Benestad}
\examples{
\dontrun{
x <- retrieve.nc("T2M_p.nc")
a <- mapField(x)
satellite(a)
}
}
\keyword{hplot}
