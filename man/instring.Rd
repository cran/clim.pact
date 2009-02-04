\name{instring}
\alias{instring}
\title{instring}
\description{Finds the position of a character in a string (character
  vector). Similar to regexpr(), but a test with regexpr() failed with
  some characters. instring() returns all position with a character
  match, whereas regexpr() only returns the first position for a pattern
  match.

  New version: in addition to finding position of one character, it also finds
  the beginning position of a given test pattern.
}
\usage{instring(c,target,case.match=TRUE)
}
\arguments{
    \item{c}{Character to look for.}
    \item{target}{string to search.}
    \item{case.match}{FALSE -> not case sensitive.}
  }
\value{vector of integers.}
\author{R.E. Benestad}
\examples{
instring("e","efile.dat")
# 1 5
regexpr("e","efile.dat")
#[1] 1
#attr(,"match.length")
#[1] 1 
# Case when regexpr() doesn't give the desired result:
regexpr(".","file.name")
#[1] 1
#attr(,"match.length")
#[1] 1
instring(".","file.name")
#[1] 5
instring("where","guess where the word where is")
# 
}
\keyword{character}
