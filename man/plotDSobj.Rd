\name{plotDSobj}
\alias{climate diagnostics}
\alias{plotDSobj}
\title{Plotting routine for objective downscaling of monthly means}
\description{Plotting routine for \code{\link{objDS}}.
}
\usage{plotDSobj(result,outdir="output",figs=c(1,2,3,4)))}

\arguments{
    \item{result}{A 'objDS' object from \code{\link{objDS}}.}
    \item{outdir}{Directory for storing EPS-files with figures}}.}
    \item{figs}{Figs to plot: 1 - time series, 2 - time series of
    residuals, 3 - distribution of residuals, 4 - rates of change.}.}
}

\value{}
}
\author{R.E. Benestad}}
\examples{
\dontrun{
library(clim.pact)
source("clim.pact/R/ds.R")
source("clim.pact/R/objDS.R")
source("clim.pact/R/catFields.R")
source("clim.pact/R/eof.R")
oslo<-getnordklim("Oslo-Blindern")
t2m <- retrieve.nc("~/data/analysis/DNMI_t2m.nc")
t2m.gcm <- retrieve.nc("~/data/mpi/mpi-gsdio_t2m.nc",
           x.rng=range(t2m$lon),y.rng=range(t2m$lat))
ds <- objDS(field.obs=t2m,field.gcm=t2m.gcm,station=oslo,lsave=FALSE,plot=FALSE)
plotDSobj(result)
}
}
\keyword{models}
\keyword{multivariate}
\keyword{ts}
\keyword{spatial}
