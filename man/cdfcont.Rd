\name{cdfcont}
\alias{cdfcont}
\title{netCDF names and dimensions.}
\description{Obtains the variable names and dimension names from a
  netCDF file. This routine provides some of the information returned by
  \code{\link{names.netCDF}}, however this function does not depend on
  the \code{netCDF} library.}
\usage{cdfcont(filename)}
\arguments{
    \item{filename}{name of netCDF file.}
 }
\value{a list object containing two character vectors: \code{vars} and
  \code{dims} holding the names of the variables and
  dimensions. \code{dims} gives the dimensions for each of the variable,
  and may contain several names within the same string, separated by a
  comma, eg. "LONS","LATS",and "LATS, LONS" (\code{strsplit("TIME, LATS,
  LONS",", ")} can be used to split up these variables into seperate names)}
\author{R.E. Benestad}}
\examples{
}
\keyword{file}
