\name{plotEOF}
\alias{plotEOF}
\title{Plot EOFs}
\description{Plots the results from an (mixed-common) EOF analysis (\code{\link{EOF}}).}
\usage{plotEOF(x,i.eof=1,nlevs=5,
                   col=c("red","blue","darkgreen","steelblue"))}
\arguments{
    \item{x}{An \code{\link{EOF}} object.}
    \item{i.eof}{EOF to plot.}
    \item{nlevs}{Contour levels.}
    \item{col}{Colour.}
    }
\value{}
\author{R.E. Benestad}}
\examples{
\dontrun{
data(eof.mc)
plotEOF(eof.mc)
x11()
data(eof.dmc)
plotEOF(eof.dmc)
}
}
\keyword{hplot}
