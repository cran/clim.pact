\name{lower.case}
\alias{lower.case}
\title{convert to lower case}
\description{Converts characters to lower case.}
\usage{lower.case(u.case)}
\arguments{\item{u.case}{Strings or arrays of strings}.}
\value{converted strings or arrays of strings.}
\author{R.E. Benestad}
\examples{
print(upper.case(c("qwerty  e","asdf  rT")))     # "QWERTY" "ASDF"
print(lower.case(c("QWERTY","ASDF")))            # "qwErty" "asdf"
print(strip(c("Hello there!","Oslo"," ","NA "))) # "Hello" "Oslo"  " "  "NA"
}
\keyword{character}
