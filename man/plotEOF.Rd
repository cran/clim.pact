\name{plotEOF}
\alias{plotEOF}
\title{Plot EOFs}
\description{Plots the results from an (mixed-common) EOF analysis
  (\code{\link{EOF}}). This function produces 3 plots: i) EOF map, ii)
  variances, and iii) PC time series. Hard copies are also produced in
  EPS format under the current working directory: "plotEOF\_1.eps",
  "plotEOF\_2.eps" and "plotEOF\_3.eps" respectively.
}
\usage{plotEOF(x,i.eof=1,nlevs=5,sub=NULL,
                   col=c("red","blue","darkgreen","steelblue"))}
\arguments{
    \item{x}{An \code{\link{EOF}} object.}
    \item{i.eof}{EOF to plot.}
    \item{nlevs}{Contour levels.}
    \item{sub}{subtitle.}
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
