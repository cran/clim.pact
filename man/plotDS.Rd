\name{plotDS}
\alias{plotDS}
\title{Plot downscaled results}
\description{Plots the results from \code{\link{DS}}.}
\usage{plotDS(ds.obj,leps=FALSE,plot.map=TRUE,plot.res=FALSE,
plot.rate=FALSE,add=FALSE,col="darkred",lwd=2,lty=1,
                   direc="output/")}
\arguments{
    \item{ds.obj}{A \code{\link{DS}} object.} 
    \item{leps}{'TRUE' produces EPS files.}
    \item{plot.map}{'TRUE' shows the spatial predictor pattern.} 
    \item{plot.res}{'TRUE' shows statistics for residuals.} 
    \item{plot.rate}{'TRUE' shows analysis of rate-of-change.}
    \item{add}{TRUE: adds a scenario to old time series plot.}
    \item{col}{Colour of scenario time series.}
    \item{lwd}{Line width of scenario time series.}
    \item{lty}{Line style of scenario time series.}
    \item{direc}{Directory for graphical output.}
    }
\value{\code{plotDS} returns a map object for the large-scale synoptic
  (anomaly) pattern, which can be plotted on a sphere
  (\code{\link{satellite}}) or coloured map (\code{\link{map}}).}
\author{R.E. Benestad}}
\examples{
data(helsinki.t2m)
data(eof.mc)
ds.helsinki<-DS(dat=helsinki.t2m,preds=eof.mc,plot=FALSE)
plotDS(ds.helsinki,leps=TRUE)
}
\keyword{hplot}