\name{julday}
\alias{julday}
\title{Converts from month, day, and year to Julian days}
\description{The function computes Julian
  days from month, day, and year. The code is based on the algortithm
  from Press et al. (1989), "Numerical Recipes  
  in Pascal", Cambridge, p. 10. See also \code{chron} and \code{date}
for similar functions. This function was included to avoid the
dependency to the \code{chron} and \code{date} packages. See also \code{\link{caldat}}.}
\usage{julday(mm,id,iyyy)}
\arguments{
    \item{mm}{month.}
    \item{id}{day.}
    \item{iyyy}{year.}
 }
\value{real}
\author{R.E. Benestad}}
\examples{
julday(1,1,1)                       # 1721424
julday(1,1,1970)                    # 2440588
julday(9,4,2003)                    # 2452887
julday(9,4,2003)-julday(1,1,1970)   # 12299
julday(9,4,2003)-julday(1,1,2003)   # 246
julday(1,1,2003)-julday(1,1,2002)   # 365
julday(1,1,2001)-julday(1,1,2000)   # 366
}
\keyword{manip}
