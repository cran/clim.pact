\name{catFields}
\alias{catFields}
\alias{fastRegrid}
\alias{testfastRegrid}
\title{catFields}
\description{Concatinates fields two different gridded sets of
  observation. The two fields must be stored on the same spatial grid,
  and the routine performs a bilinear spatial interpolation to place the
  data on the same grid. Observations/data for representing values at n
  different locations at a given time (t) can be described in terms of a
  vector
  \deqn{\vec{x}(t) = [x_1, x_2, ... x_n].}{x(t)=[x_1, x_2,
    ... x_n].}
  The data set consists of a time series of vectors which can be
  represented by the means of matrices
  \deqn{X = [\vec{x}(t_1), \vec{x}(t_2), ... \vec{x}(t_n)].}{
    X=[x(t_1), x(t_2), ... x(t_n)].}
  Two different sets of observations can be represented by two matrices
  Y and Z with dimensions k x n and k x m respectively (k is the number
  of spatial points, whereas n and mIindicate the number of observations
  in time). The information in these two data sets are combined
  combining  the two matrices using rbind. The major difference between
  this routine and rbind is that this routine takes care of all the
  'house keeping' in terms of grid, time and variable information. 

  cat.field can be used to process single fields by setting
  'field.2=NULL'. This option allows for interpolation and extraction of
  sub-regions or sub-intervals, removing the mean values, and selecting
  a paricular month or season.

  The output from cat.fields can by further analysed in  \code{\link{EOF}}. By
  using a concatination of two fields of similar data, eg observed and
  simulated sea level pressure (SLP), it is possible to carry out a
  common EOF analysis. The application of \code{\link{DS}} to the EOFs of
  concatinated fields provides an analysis similar to the common EOH
  method described in Benestad (2001), "A comparison between two
  empirical downscaling strategies", \emph{Int. J. Climatology}, \bold{vol 21},
  1645-1668, DOI 10.1002/joc.703.

  \code{fastRegrid} is a quick way to re-grid field object, and uses
  \code{\link{EOF}} to decompose the field into a small number of
  spatial patterns (EOFs), which are then re-gridded, and subsequently
  transformed back to a field object by calling \code{\link{EOF2field}}.
}  
\usage{
catFields(field.1,field.2=NULL,lat=NULL,lon=NULL,plot.interp=FALSE,
           interval.1=NULL,interval.2=NULL,mon=NULL,demean=TRUE,silent=FALSE,
                      fastregrid=FALSE,neofs=20)
fastRegrid(field.1,lat=NULL,lon=NULL,mon=NULL,neofs=20,
                      silent=TRUE,plot.interp=FALSE)
testfastRegrid()
}
\arguments{
    \item{field.1}{A 'field.object'. }
    \item{field.2}{A 'field.object'. A 'field.2=NULL' processes single fields.}
    \item{lat}{Latitudes to extract. If NULL, use the latitudes from
      the first field. Otherwise interpolate both fields to latitudes.  }
    \item{lon}{Longitudes to extract. See 'lat'. }
    \item{plot.interp}{Flag: 'TRUE' plots the interploation results - 
    Used for checking interpolation. }
    \item{interval.1}{Extract the time interval for the 1st field. }
    \item{interval.2}{Extract the time interval for the 2nd field. }
    \item{mon}{Calendar month or season to extract. eg January or DJF. }
    \item{demean}{Flag: 'TRUE' subtracts the mean values. This flag
      should be set to 'FALSE' if for instance two time slices are
      concatinated and the object is to investigate the mean change
      between these periods (see examples in \code{\link{eof.dc}} or
      \code{\link{eof.dmc}}).}
    \item{silent}{FALSE: verbose output.}
    \item{fastregrid}{Option in catFields - use \code{fastRegrid}}
    \item{neofs}{Number of EOF-modes to retain, determining the degree
  of details to be retained.}
  }
\value{A 'field.object'.}
\author{R.E. Benestad}
\examples{
\dontrun{
library(clim.pact)
x.1 <- retrieve.nc("/home/kareb/data/ncep/ncep_t2m.nc",
                   x.rng=c(-60,40),y.rng=c(50,75))
x.2 <- retrieve.nc("/home/kareb/data/ncep/ncep_slp.nc",
                   x.rng=c(-60,40),y.rng=c(50,75))
print(x.1$v.name)

print("Read GCM predictor data.")
X.1 <- retrieve.nc("data/mpi-gsdio_t2m.nc",
                   x.rng=c(-60,40),y.rng=c(50,75))
X.2 <- retrieve.nc("data/mpi-gsdio_slp.nc",
                   x.rng=c(-60,40),y.rng=c(50,75))
print(X.1$v.name)
print("Cat fields.")
xX.1 <- catFields(x.1,X.1,interval.1=c(1958,1998),interval.2=c(1958,2050))
xX.2 <- catFields(x.2,X.2,interval.1=c(1958,1998),interval.2=c(1958,2050))
xX <- mixFields(xX.1,xX.2,mon=1,
                 interval=c(1900,2050))
print("EOF")
eof.c <- eof(xX.1,mon=1)
eof.mc <- eof(xX,mon=1)

testfastRegrid()

data(DNMI.slp)
slp <- fastRegrid(DNMI.slp,lat=seq(30,80,by=1),lon=seq(-90,60,by=1),mon=1)
}
}
\keyword{manip}
\keyword{ts}
