\name{satellite}
\alias{polar stereographic}
\alias{satellite}
\title{Satellite view / polar stereographic}
\description{Produces polar stereographic maps / satellite views.}
\usage{satellite(map.obj,col="black",lwd=2,lty=1,add=FALSE,
                      las = 1,lon.0=NULL,lat.0=NULL,method="normal",
                      ni=100,nj=100, n.nearest=4,max.dist=3,landdata=addland2)}
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
