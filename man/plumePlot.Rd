\name{plumePlot}
\alias{plumePlot}
\title{Plot downscaled time series as plumes}
\description{Retreives an ensemble of time series from \code{\link{DS}}
  and plots these as a plume. The routine retrieves the ds objects from disc.}
\usage{plumePlot(ds.name.list=NULL,location,mon,direc="output",
                         t.rng=c(1850,2074),r2.th=50,p.th=0.05,
                         col="darkred",lwd=2,lty=1)}
\arguments{
  \item{ds.name.list}{A list of file names holding the ds objects.}
  \item{location}{Name of location to plot.}
  \item{mon}{Month to plot.}
  \item{direc}{The subdirectory in which the ds objects are stored.}
  \item{t.rng}{Time interval to plot.}
  \item{r2.th}{R-squared theshold: only use scenarios that account for
      equal to or more than \code{r2.th} of the variance in \%.}
  \item{p.th}{p-value threshold: only use scenarios that have p-values
      equal to or less than \code{p.th}.}
  \item{col}{Colour of scenario time series.}
  \item{lwd}{Line width of scenario time series.}
  \item{lty}{Line style of scenario time series.}
 }
\value{}
\author{R.E. Benestad}}
\examples{
\dontrun{
ds.list<-avail.ds()
plumePlot(ds.list,location="OSLO-BLINDERN",mon=1)
}
}
\keyword{hplot}