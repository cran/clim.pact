\name{objDS}
\alias{climate analysis}
\alias{objDS}
\alias{downscaling}
\title{Objective downscaling of monthly means}
\description{'Objective' downscaling based on \code{\link{DS}}. The
function selects region according to a correlation analysis, setting
the borders where the correlation with the station series become
zero. A fit to a truncated Fourier expansion is used to describe
profiles of correlation coefficients in zonal and meridional
directions from where the station is located.

Version 2.1-5: smallest spatial domain for objDS is set to +- 10deg N 
and E of the station location in order to avoid problems with dimensions 
if one only contains one data point (the domain becomes a line instead
of a rectangle).
}
\usage{objDS(field.obs,field.gcm,station,plot=TRUE,positive=NULL,
          mon=NULL,direc="output/",cal.id=NULL,
          ldetrnd=TRUE,i.eofs=seq(1,8,by=1),ex.tag="",
          method="lm",leps=FALSE,param="t2m",
          plot.res=FALSE,plot.rate=FALSE,xtr.args="",
          swsm="step",predm="predict",lsave=FALSE,rmac=TRUE,
          silent=FALSE,qualitycontrol=TRUE)}

\arguments{
    \item{field.obs}{ The gridded observation predictor \code{\link{retrieve.nc}}.}
    \item{field.gcm}{ The climate simulation predictor \code{\link{retrieve.nc}}.}
    \item{station}{ A climate.station object (\code{\link{station.obj}} or 
        \code{\link{station.obj.dm}}). 
        [e.g. from \code{getnacd}, \code{getnordklim} or \code{station.obj}].}
    \item{plot}{'TRUE' produces figures.}
    \item{positive}{'TRUE': only consider the region where correlations
      are positive (important for temperature predictors).}
    \item{mon}{ month or season to downscale,loops though  the 12
      calendar months if NULL.}
    \item{direc}{ name of directory inwhich the output is dumped
          (e.g. figures, tables).}
    \item{cal.id}{ID tag used for calibration. By default use the 
	first field (\code{\link{catFields}}) for calibration.}
    \item{ldetrnd}{ F for no detrending; T for removing linear trends
	before model calibration.}
    \item{i.eofs}{ select which EOFs to include in the setp-wise screening.}
    \item{ex.tag}{Extra labelling tag for file names for experiments.}
    \item{method}{Sets the method to use for regression. Method is set
      to "lm" by default, but "anm" allows the incorporation of an
      analog model, see \code{\link[anm]{anm}}. "anm.weight" weights the
    principal components according to the eigenvalues, whereas "anm"
    uses unweighted series.}
    \item{leps}{'TRUE' produces EPS figures (files).}
    \item{param}{Name of parameter (for plot labels).}
    \item{plot.res}{'TRUE' shows statistics for residuals.}
    \item{plot.rate}{'TRUE' shows analysis of rate-of-change.}
    \item{xtr.args}{Extra/additional arguments in the formula.}
    \item{swsm}{Step-wise screening method, default=='step'; 'none'
      skips stepwise sceeening.}
    \item{predm}{Prediction method, default is "predict"}
    \item{lsave}{TRUE -> saves the result on disc}.
    \item{rmac}{TRUE -> subtracts (removes) the annual cycle in station
      data.}
    \item{silent}{TRUE -> no output to screen.}
    \item{qualitycontrol}{TRUE: perform a qualitry control consisting of comparing the smoothness of the monthly trend estimates throughout the year and repeat the downscaling with a smaller domain if adjacent trend estimates vary significantly (diff greater than 3*variance of rates).}
}

\value{An 'objDS' object -  a list of objects:
  \tabular{rl}{
   station \tab a 'station' object  (see \code{\link{getnacd}}) \cr
   Jan    \tab a 'ds' (see \code{\link{DS}}) object \cr
   Feb    \tab a 'ds' object \cr
   Mar    \tab a 'ds' object \cr
   Apr    \tab a 'ds' object \cr
   May    \tab a 'ds' object \cr
   Jun    \tab a 'ds' object \cr
   Jul    \tab a 'ds' object \cr
   Aug    \tab a 'ds' object \cr
   Sep    \tab a 'ds' object \cr
   Oct    \tab a 'ds' object \cr
   Nov    \tab a 'ds' object \cr 
   Dec    \tab a 'ds' object \cr 
  }
}
\author{R.E. Benestad}}
\examples{
\dontrun{
library(clim.pact)
oslo<-getnordklim("Oslo-Blindern")
slp.obs <- retrieve.nc("ncep_slp.nc")         # Get gridded observations/analysis from NCEP
slp.gcm <- retrieve.nc("EH4OPYC_B2_slp.nc")   # Get results from climate models
ds <- objDS(field.obs=slp.obs,field.gcm=slp.gcm,station=oslo)

t2m <- retrieve.nc("~/data/analysis/DNMI_t2m.nc")
t2m.gcm <- retrieve.nc("~/data/mpi/mpi-gsdio_t2m.nc",
           x.rng=range(t2m$lon),y.rng=range(t2m$lat))
ds <- objDS(field.obs=t2m,field.gcm=t2m.gcm,station=oslo)
}
}
\keyword{models}
\keyword{multivariate}
\keyword{ts}
\keyword{spatial}
