\name{strip}
\alias{strip}
\title{String operation functions}
\description{The function strips off trailing space (strips the
  strings by cutting off at the first space).}
\usage{strip(in.str)}
\arguments{\item{in.str}{Strings or arrays of strings.}}
\value{converted strings or arrays of strings.}
\author{R.E. Benestad}}
\examples{
       print(upper.case(c("qwerty  e","asdf  rT")))       #  "QWERTY" "ASDF"
       print(lower.case(c("QWERTY","ASDF")))              #  "qwErty" "asdf"
       print(strip(c("Hello there!","Oslo"," ","NA ")))   # "Hello" "Oslo"  " "     "NA"
}
\keyword{character}
