\name{retrieve.nc}
\alias{retrieve.nc}
\alias{field.object}
\alias{fixField}
\alias{monthly}
\alias{test.retrieve.nc}
\title{Retrieve data from a netCDF file}
\description{\code{retrieve.nc} reads a netCDF file and picks out vectors that look like
  lngitude, latitude and time. Returns the first 3-D field in the file.
  See also \code{\link{cdfextract}} (large data files). 

This version uses \code{\link{cdfcont}} to obtain vital meta data for
handling the data in the netCDF file and constructing a 'field' object.

This routine reequres the \code{\link[ncdf]{ncdf}}-package is currently only available for Linux systems at the CRAN web sites, but the Windows version can can be obtained from \url{http://www.stats.ox.ac.uk/pub/RWin/}.

\code{fixField} is a function for eg changing the time stamp if the given time origin is wrong.
}

\usage{
retrieve.nc(filename=file.path("data","air.mon.mean.nc"),v.nam="AUTO",
            l.scale=TRUE,greenwich=TRUE,silent=FALSE,
            x.nam="lon",y.nam="lat",z.nam="lev",t.nam="tim",
            x.rng=NULL,y.rng=NULL,z.rng=NULL,t.rng=NULL,
            force.chron=TRUE,force365.25=FALSE,regular=TRUE,
            daysayear=365.25,forceBC=TRUE,
            use.cdfcont=FALSE,torg=NULL,t.unit=NULL,miss2na=TRUE)
fixField(x,torg=NULL,t.unit=NULL,scal=NULL,offs=NULL, 
         x.rng=NULL,y.rng=NULL,z.rng=NULL,t.rng=NULL,greenwich=TRUE)
test.retrieve.nc(filename="sst.wkmean.1981-1989.nc")
monthly(x,param="t2m",method="mean")
}

\arguments{
    \item{filename}{name of netCDF file.}
    \item{v.nam}{name of variable. "AUTO" -> smart search.}
    \item{l.scale}{'TRUE' uses scaling.factor and add.offset.}
    \item{greenwich}{'TRUE' centres maps on Greenwhich meridian (0 deg E).}
    \item{x.nam}{Name of x-dimension.}
    \item{y.nam}{Name of y-dimension.}
    \item{z.nam}{Name of z-dimension.}
    \item{t.nam}{Name of time-axis.}
    \item{x.rng}{Region to extract.}
    \item{y.rng}{Region to extract.}
    \item{z.rng}{Region to extract.}
    \item{t.rng}{Time interval to extract.Numerical values are used to identify indeces, e.g. \code{as.numeric(1)} refers to first field, \code{as.numeric(2)} the second field, etc. Character strings, on the other hand, refers to date. E.g. "1-Jan-1990", or "1990" (see \code{\link{datestr2num}} for various formats).}
    \item{force.chron}{Check for monotonic chronological order (no
      jumping back and forth in time).}
    \item{force365.25}{TRUE forces a natural 365.25 day year as opposed to a 360-day model year. '-1' forces a 360-day year (commonly used in climate modelling).}
    \item{regular}{TRUE for regular spacing in time (.i.e. no skipping, but one field every month or one field every day).}
    \item{daysayear}{Number of days in the year on average.}
    \item{forceBC}{TRUE for not accepting year 0 (see e.g. Press et al. (1989), Numerical Recepies in Pascal, Cambridge).}
    \item{x}{a field object.}
    \item{scal}{scaling factor}
    \item{offs}{constant offset}
    \item{silent}{To turn off printing to console.}
    \item{param}{Parameter to convert to monthly value}
    \item{use.cdfcont}{Flag for Linux versions only: if TRUE use old lines calling 'cdfcont()'}
    \item{torg}{Time origin, such as the 'time\_origin' attribute in netCDF
      files. e.g. '15-Dec-1949'. A NULL value (default) will try to detect from
      the file header.}
    \item{t.unit}{Time unit, similar to the 'time\_unit' attribute in
      netCDF files. e.g. 'day'. A NULL value (default) will try to detect from
      the file header.}
    \item{method}{Function to use when computing monthly values
      (e.g. mean, sd, max, min)}
    \item{miss2na}{FLAG: TRUE - set missing value from netCDF meta data to NA}
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
\author{R.E. Benestad}
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

# For reading the IPCC FoAR netCDF files that uses a 365-day year (no leap years)
#  and starts on time count year 0:
gcm <- retrieve.nc(fname,v.nam="tas",x.rng=c(-50,50),y.rng=c(30,75),forceBC=FALSE)
}
}
\keyword{file}
