\name{r2cdf}
\alias{r2cdf}
\title{Save as netCDF file.}
\description{Saves a field object, map or eof object as a netCDF file by
  creating a CDF file (ASCII) and then calling 'ncgen' through \code{system()}.}
\usage{r2cdf(filename,x,missing=-999.99,cleanup=TRUE,
                  ofs=NULL,scal=NULL)}
\arguments{
  \item{filename}{Name of file.}
  \item{x}{a Field, map, or eof object}
  \item{missing}{Missing value}
  \item{cleanup}{TRUE removes the CDF-file after the call.}
  \item{ofs}{offset: 'NULL' -> automatic setting}
  \item{scal}{scaling factor: 'NULL' -> automatic setting}
 }
\value{}
\author{R.E. Benestad}
\examples{
\dontrun{
# Save EOFs as netCDF (use ncview or Ferret to view)
data(eof.slp)
r2cdf("test.nc",eof.slp)
#
#
slp <- retrieve.nc("data/DNMI_slp.nc")
mslp <- meanField(slp)
r2cdf("test.nc",mslp)
r2cdf("test.nc",slp)
#
slp <- cdfextract("data/nmc_slp.nc","slp",x.rng=c(-80,40),y.rng=c(20,75),
                  t.rng=c(times[is],times[is]+499),plot=FALSE)
r2cdf("test.nc",slp)
#
data(oslo.t2m)
map <- composite.field(slp,oslo.t2m)
r2cdf("test.nc",map)
#
Xdum=list(dat=slp$dat[1:10,,],lon=slp$lon,lat=slp$lat,tim=slp$tim[1:10],
          lev=NULL,v.name=slp$v.name,attributes=slp$attributes)
class(Xdum)="field"
r2cdf("test.nc",Xdum)
#
# The definition of a 'field' object is:
ny<-length(slp$lat); nx<-length(slp$lon)
slp  <- list(dat=slp$dat,lon=slp$lon,lat=slp$lat,tim=slp$tim,lev=slp$lev,
             v.name=slp$v.nam,id.x=slp$id.x,id.t=slp$id.t,
             yy=slp$yy,mm=slp$mm,dd=slp$dd,n.fld=1,
             id.lon=rep(slp$v.name,nx),id.lat=rep(slp$v.name,ny),
             attributes=dat.att)
 class(slp) <- c("field")

}
}
\keyword{file}
