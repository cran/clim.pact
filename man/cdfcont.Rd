\name{cdfcont}
\alias{cdfcont}
\title{netCDF names and dimensions.}
\description{Obtains the variable names and dimension names from a
  netCDF file. This routine provides some of the information returned by
  the old \code{netCDF}-package function  \code{names.netCDF},
  however this function does not depend on
  the \code{netCDF} package. If it doesn't work, make sure to set the
  path to include the location of \code{ncdump} and \code{ncgen} (two
  netCDF executables, usually under a netcdf direcory). Alternatively,
  set the \code{path}.}
\usage{cdfcont(filename,path="",method=NULL)}
\arguments{
    \item{filename}{name of netCDF file.}
    \item{path}{Path to the location where \code{ncdump} and \code{ncgen} are located.}
    \item{method}{The name of the function making calls to the
      system. In Linux, it is 'system', in Windows 'shell'. If set to
      NULL, the code tries to decide between these two by probing using \code{\link{Sys.info}}}
 }
\value{a list object containing two character vectors: \code{vars} and
  \code{dims} holding the names of the variables and
  dimensions. \code{dims} gives the dimensions for each of the variable,
  and may contain several names within the same string, separated by a
  comma, eg. "LONS","LATS",and "LATS, LONS" (\code{strsplit("TIME, LATS,
  LONS",", ")} can be used to split up these variables into seperate
  names).
  'time.origin' returns the time\_origin, 'add.offset' gives the 'add\_offset'
  attribute, and 'scale.factor' gives the scale\_factor.}

\author{R.E. Benestad}
\examples{
}
\keyword{file}
