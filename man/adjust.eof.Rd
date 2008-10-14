\name{adjustEOF}
\alias{adjustEOF}
\title{Adjust common EOFs}
\description{Forces same std and mean in the subsegments of common PCs.} 
\usage{adjustEOF(x)}
\arguments{
    \item{x}{a 'eof' object (see \code{\link{EOF}}).}
 
  }
\value{an 'eof' object}
\author{R.E. Benestad}
\examples{
data(eof.c)
eof <- adjustEOF(eof.c)
plotEOF(eof)
}
\keyword{manip}
