\name{num2str}
\alias{num2str}
\alias{strings}
\alias{characters}
\title{Convert numbers to string and format}
\description{Convert numbers to string and format. Similar to FORTRAN
  'F8.2' format statement.}
\usage{num2str(x,dec=2,f.width=NULL,d.point=".")}
\arguments{ \item{x}{Real numbers.}
            \item{dec}{number of decmal points.}
            \item{f.width}{width of field.}
            \item{d.point}{character marking the decimal point.}
\value{string.}
\author{R.E. Benestad}}
\examples{
print(num2str(c(1,23.4282,-3.14),dec=3))
#[1] "1.000"  "23.428" "-3.140"
}
\keyword{character}