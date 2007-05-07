\name{datestr2num}
\alias{datestr2num}
\title{datestr2num}
\description{Converts a date string to number.}
\usage{datestr2num(datestr,vec=TRUE)}
\arguments{
    \item{datestr}{A date string (character). Types of stings may be "01-Jan-1980","1-1-1980", "01-01-1980", "1980-1-1", "1980-Jan-1"}
    \item{vec}{TRUE returns a vector \code{c(year,month,day)}; FALSE returns \code{year + month/12 + day/31}.}
 }
\value{A number or a vector}
\author{R.E. Benestad}
\examples{
datestr2num("01-Jan-1980")
datestr2num("1-1-1980")
datestr2num("01-01-1980")
datestr2num("1980-1-1")
datestr2num("1980-Jan-1")
}
\keyword{manip}
