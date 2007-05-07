\name{reverse}
\alias{reverse}
\alias{reverse.sort}
\title{Reverse}
\description{Reverses the order of a vector.\code{reverse.sort} returns
  a sorted vector in reverse order. }
\usage{reverse(x)
reverse.sort(x)
}
\arguments{
    \item{x}{a vector}
  }
\value{a vector.}
\author{R.E. Benestad}
\examples{
reverse(c(1,3,5,7,2,4,6,8))      # 8 6 4 2 7 5 3 1
reverse.sort(c(1,3,5,7,2,4,6,8)) # 8 7 6 5 4 3 2 1
}
\keyword{manip}
