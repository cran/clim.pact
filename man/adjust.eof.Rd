\name{adjust.eof}
\alias{adjust.eof}
\title{Adjust common EOFs}
\description{Forces same std and mean in the subsegments of common PCs.} 
\usage{adjust.eof(x)}
\arguments{
    \item{x}{a 'eof' object (see \code{\link{EOF}}).}
 
  }
\value{an 'eof' object}
\author{R.E. Benestad}}
\examples{
data(eof.c)
eof <- adjust.eof(eof.c)
plotEOF(eof)
}
\keyword{manip}
