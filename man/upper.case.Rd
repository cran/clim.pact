\name{upper.case}
\alias{upper.case}
\title{convert to UPPER CASE}
\description{Converts characters to UPPER CASE.}
\usage{upper.case(u.case)}
\arguments{\item{u.case}{Strings or arrays of strings.}}
\value{converted strings or arrays of strings.}
\author{R.E. Benestad}
\examples{
print(upper.case(c("qwerty  e","asdf  rT")))     #  "QWERTY" "ASDF"
print(lower.case(c("QWERTY","ASDF")))            #  "qwErty" "asdf"
print(strip(c("Hello there!","Oslo"," ","NA "))) # "Hello" "Oslo" " "  "NA"
}
\keyword{character}
