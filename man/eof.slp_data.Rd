\name{eof.slp}
\alias{eof.slp}
\title{EOF of NCEP reanalysis SLP.}
\description{EOF of NCEP reanalysis January mean sea level pressure
}
\usage{data(eof.slp)}
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
The common EOF was produced using  \code{\link{EOF}}, with the combined January
2-meter air temperature data field from National Center for
Environmental Prediction (NCEP; USA) reanalysis (Kalnay et al., (1996) "The NCEP/NCAR 40-Year Reanalysis Project",  \emph{Bul. Am. Met. Soc.}, \bold{vol 77}, no 3, 437-471; e.g. see URL: \url{http://www.cdc.noaa.gov/index.html}).
}
\references{Reference to methodology: R.E. Benestad (2001), "A comparison
  between two empirical downscaling strategies",
  \emph{Int. J. Climatology}, \bold{vol 210}, pp.1645-1668. [DOI
  10.1002/joc.703].
  }
\examples{
library(clim.pact)
data(eof.slp)
}
\keyword{datasets}
