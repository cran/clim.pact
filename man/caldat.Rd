\name{caldat}
\alias{caldat}
\title{Converts Julian days to month, day, and year}
\description{The function computes month, day, and year from Julian
  days. The code is based on the algortithm from Press et al. (1989), "Numerical Recipes 
  in Pascal", Cambridge, p. 13. See also \code{chron} and \code{date}
for similar functions. This function was included to avoid the
dependency to the \code{chron} and \code{date} packages. See also \code{\link{julday}}.}
\usage{caldat(julian)}
\arguments{
    \item{julian}{The Julian day from 1-1-1.}
 }
\value{a list with: month, day, and year}
\author{R.E. Benestad}
\examples{
caldat(1)        # month=1, day=2, year=-4713
caldat(1721424)  #       1,     1,          1
caldat(2440588)  #       1,     1,       1970
caldat(2452887)  #       9,     4,       2003
}
\keyword{manip}
