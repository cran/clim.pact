\name{addClim}
\alias{addClim.station}
\alias{preClim.station}
\title{Add climatology}
\description{Adds the climatological values to anomalies. The inverse of \code{\link{anomaly.station}}.

}
\usage{
addClim.station(x,param=c("t2m","precip"))
preClim.station(x,dd=NULL,mm=NULL,yy=NULL,param=1)
}
\arguments{
    \item{x}{A station object or a vector.}
    \item{mm}{month}
    \item{dd}{day.}
    \item{yy}{year.}
    \item{param}{For daily objects: identifies the parameters (addClim.station) or determines which model to use (preClim.station).}
 }
\value{A station object (addClim.station) or a vector (preClim.station).}
\author{R.E. Benestad}
\examples{
data(oslo.dm, envir=environment())
oslo.dma <- anomaly.station(oslo.dm)
oslo.dmac <- addClim.station(oslo.dma)
plot(oslo.dm$yy+(oslo.dm$mm-1)/12+oslo.dm$dd/365.25,oslo.dm$t2m,
     type="l",lwd=3,col="grey",xlim=c(1980,1985))
lines(oslo.dm$yy+(oslo.dm$mm-1)/12+oslo.dm$dd/365.25,
      oslo.dmac$t2m,lty=2)
lines(oslo.dm$yy+(oslo.dm$mm-1)/12+oslo.dm$dd/365.25,
      preClim.station(oslo.dma, oslo.dm$dd, oslo.dm$mm,
      oslo.dm$yy),col="red",lwd=2)
}
\keyword{manip}
