\name{addland}
\alias{addland}
\title{Add land contours to map.}
\description{The function superimposes land contours on a map.}
\usage{addland(col="grey50",lwd=1) }
\arguments{
 \item{col}{colour}
 \item{lwd}{line width}
}
\author{R.E. Benestad}}
\examples{
plot(c(-90,90),c(0,80),type="n")
addland()
grid()
}
\keyword{hplot}
