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
\author{R.E. Benestad}}
\examples{
}
\keyword{file}
