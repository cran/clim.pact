\name{mod}
\alias{mod}
\title{Modulus of a division.}
\description{Returns the modulus of a division: returns the remainder of
the expression x/y.}
\usage{mod(x,y)}
\arguments{
  \item{x}{nominator (integer).}
  \item{y}{denominator (integer).}
}
\value{integer value from 0 .. (y-1)}
\author{R.E. Benestad}}
\examples{
mod(101,10)   # 1
mod(4,12)     # 4
mod(123,12)   # 3
}
\keyword{arith}

