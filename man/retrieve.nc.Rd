\name{retrieve.nc}
\alias{retrieve.nc}
\alias{field.object}
\title{Retrieve data from a netCDF file}
\description{Reads a netCDF file and picks out vectors that look like
  lngitude, latitude and time. Returns the first 3-D field in the file.
  See also \code{\link{cdfextract}} (large data files).}
\usage{retrieve.nc(f.name="data/ncep_t2m.nc",v.nam="AUTO",
                        l.scale=TRUE,greenwich=TRUE,
                        x.nam="lon",y.nam="lat",z.nam="lev",t.nam="tim",
                        x.rng=NULL,y.rng=NULL,t.rng=NULL)
}
\arguments{
    \item{f.name}{name of netCDF file.}
    \item{v.nam}{name of variable. "AUTO" -> smart search.}
    \item{l.scale}{'TRUE' uses scaling.factor and add.offset.}
    \item{greenwich}{'TRUE' centres maps on Greenwhich meridian (0 deg E).}
    \item{x.nam}{Name of x-dimension.}
    \item{y.nam}{Name of y-dimension.}
    \item{z.nam}{Name of z-dimension.}
    \item{t.nam}{Name of time-axis.}
    \item{x.rng}{Region to extract.}
    \item{y.rng}{Region to extract.}
    \item{t.rng}{Time interval to extract.}
  }
  \value{A "field.object" list:
    \tabular{ll}{
      dat \tab a 3-D matrix with the data. \cr
      lon \tab a vector of longitudes.  \cr
      lat \tab a vector of latitudes.  \cr
      tim \tab a vector of times from time.0 (see attributes). \cr
      lev \tab a vector of vertical levels (NULL for single level). \cr
      v.name \tab variable name.\cr
      id.x \tab ID labels for the spatial grid (for mixed fields, see \code{\link{mixFields}}). \cr
      id.t \tab ID labels for the time axis (for combined fields). \cr
      yy \tab a vector of years corresponding to \code{tim}. \cr
      mm \tab a vector of months corresponding to \code{tim}. \cr
      dd \tab a vector of days corresponding to \code{tim}. \cr
      n.fld \tab number of fields (for mixed fields, see
      \code{\link{mixFields}}). \cr
      id.lon \tab ID labels along the longitudes (for mixed fields, see
      \code{\link{mixFields}}). \cr
      id.lat \tab ID labels along the latitudes (for mixed fields, see
      \code{\link{mixFields}}). \cr
      }.}
\author{R.E. Benestad}}
\examples{
\dontrun{
X.1 <- retrieve.nc("data/mpi-gsdio_t2m.nc",
                   x.rng=c(-60,40),y.rng=c(50,75))
X.2 <- retrieve.nc("data/mpi-gsdio_slp.nc",
                   x.rng=c(-60,40),y.rng=c(50,75))

# The definition of a 'field' object is:
ny<-length(slp$lat); nx<-length(slp$lon)
slp  <- list(dat=X.2$dat,lon=X.2$lon,lat=X.2$lat,tim=X.2$tim,lev=X.2$lev,
             v.name=X.2$v.nam,id.x=X.2$id.x,id.t=X.2$id.t,
             yy=X.2$yy,mm=X.2$mm,dd=X.2$dd,n.fld=1,
             id.lon=rep(X.2$v.name,nx),id.lat=rep(X.2$v.name,ny),
             attributes=dat.att)
 class(slp) <- c("field")

}
}
\keyword{file}
