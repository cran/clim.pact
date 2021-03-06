\name{delta}
\alias{delta}
\alias{delta function}
\title{Delta function}
\description{The dirac delta function: delta (i,j) returns 1 if i==j, 0
  otherwise.}
\usage{delta(i,j)}
\arguments{
    \item{i}{first index.}
    \item{j}{second index.}
 }
\value{0 or 1}
\author{R.E. Benestad}
\examples{
zero <- delta(1,2)
one <- delta(2,2)
}
\keyword{manip}
