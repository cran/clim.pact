\name{retrieve.nc}
\alias{retrieve.nc}
\title{Retrieve data from a netCDF file}
\description{Reads a netCDF file and picks out vectors that look like
  lngitude, latitude and time. Returns the first 3-D field in the file.}
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
      lon \tab a vector of longitudes  \cr
      lat \tab a vector of latitudes  \cr
      tim \tab a vector of times from time.0 (see attributes) \cr
      dat \tab a 3-D matrix with the data \cr
      }.}
\author{R.E. Benestad}}
\examples{
\dontrun{
X.1 <- retrieve.nc("data/mpi-gsdio_t2m.nc",
                   x.rng=c(-60,40),y.rng=c(50,75))
X.2 <- retrieve.nc("data/mpi-gsdio_slp.nc",
                   x.rng=c(-60,40),y.rng=c(50,75))
}
}
\keyword{models}
