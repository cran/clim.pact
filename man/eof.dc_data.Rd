\name{eof.dc}
\alias{eof.dc}
\alias{eof.dc_data}
\title{Daily common EOF.}
\description{Common EOFs for daily December-February 2-meter temperature
  (T(2m)). 
}
\usage{data(eof.dc)}
\format{
  \tabular{ll}{
    EOF \tab EOF patterns. \cr
    W  \tab Eigen values. \cr
    PC \tab Principal components of common PCA. \cr
    n.fld \tab Number of different predictors (see
    \link{mixFields}). \cr
    tot.var \tab Sum of all W squared. \cr
    id.t \tab Time labels for the fields (see
    \link{catFields}) - used in \link{DS}. \cr
    id.x \tab Spatial labels for the fields (see
    \link{mixFields}) - used in \link{plotEOF}. \cr
    id.lon \tab Spatial labels for the fields (see
    \link{mixFields}) - used in \link{plotEOF}. \cr
    id.lat \tab Spatial labels for the fields (see
    \link{mixFields}) - used in \link{plotEOF}. \cr
    region \tab Describes the region analysed. \cr
    tim \tab Time information (usually redundant). \cr
    lon \tab Longitudes associated with EOF patterns. \cr
    lat \tab Latitudes associated with EOF patterns. \cr
    var.eof \tab Fractional variances associated with EOF patterns. \cr
    yy \tab years. \cr
    mm \tab months. \cr
    dd \tab days. \cr
    v.name \tab Name of element. \cr
    c.mon \tab Month-season information. \cr
    f.name \tab File name of original data. \cr
  }
  }
\source{
The common EOF was produced using  \code{\link{EOF}}, with the combined
December-February (DJF) 2-meter air temperature data field from European
Centre for Medium-Range Weather Forecasts (ECMWF; UK) reanalysis ( see URL:
\url{http://www.ecmwf.int/}) and HIRHAM dynamically dowwnscaled
scenarios from the ECHAM4-GSDIO scenario (Max-Planck Institute for Meteorology, Hamburg, Germany; URL: \url{http://www.mpimet.mpg.de/}). The region is 5E - 25E, 58N - 65N.
}
\references{Reference to methodology: R.E. Benestad (2001), "A comparison
  between two empirical downscaling strategies", \emph{Int. J. Climatology}, \bold{vol 210}, pp.1645-1668. [DOI 10.1002/joc.703].}
\examples{
#The EOFs were produced using the following code:
library(clim.pact)
\dontrun{
x.1.dm<-retrieve.nc("/data1/era15/ERA-15_t2m.nc",x.rng=c(5,25),y.rng=c(58,65))
X.1.dm<-retrieve.nc("/data1/hirham/T2M_198001-199912.nc",x.rng=c(5,25),
                    y.rng=c(58,65))
Y.1.dm<-retrieve.nc("/data1/hirham/T2M_203001-204912.nc",x.rng=c(5,25),
                    y.rng=c(58,65))
Y.1.dm$yy <- Y.1.dm$yy + 50
# It is important that demean=FALSE when concatinating the two time slices
# from the model simulations, if a study of climate change is the objective.
xX.1.dm <- cat.fields(X.1.dm,Y.1.dm,demean=FALSE)
xX.1.dm <- cat.fields(x.1.dm,xX.1.dm)
eof.dc <- eof(xX.1.dm,mon=1)
}
# To read the data:
data(eof.dc)
}
\keyword{datasets}
