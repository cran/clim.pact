\name{cdfextract}
\alias{cdfextract}
\title{Extract a subfield from a netCDF file.}
\description{A slow routine that extraxcts a subfield. This routine is
  suitable for reading subsections of large data files that are too big
  for \code{\link{retrieve.nc}}.}
\usage{cdfextract(filename,varname,x.rng=NULL,y.rng=NULL,t.rng=NULL,
                  greenwich=TRUE,x.nam="lon",y.nam="lat",t.nam="tim",
                       plot=FALSE,l.scale=TRUE)}
\arguments{
    \item{filename}{name of netCDF file.}
    \item{varname}{name ov variable.}
    \item{x.rng}{X-range (in degrees East): c(min,max).}
    \item{y.rng}{Y-range (in degrees North): c(min,max).}
    \item{t.rng}{T-range (in units stored in 'tim'): c(min,max).}
    \item{greenwich}{TRUE: longitude runs from -180E to 180E.}
    \item{x.nam}{name of X-dimension.}
    \item{y.nam}{name of Y-dimension.}
    \item{t.nam}{name of T-dimension.}
    \item{plot}{TRUE for plotting.}
    \item{l.scale}{scale field by offset and scaling-factor.}
 }
\value{a field object. Also see \code{\link{retrieve.nc}}. Saves the
  extracted data in a netCDF file called "cdfextract.nc" under current
  working directory (see \code{\link{r2cdf}}).}
\author{R.E. Benestad}}
\examples{
\dontrun{
slp <- cdfextract("data/nmc_slp.nc","slp",x.rng=c(-80,40),y.rng=c(20,75),
                  t.rng=c(times[is],times[is]+499),plot=FALSE)
}
}
\keyword{file}
