\name{patternIndex}
\alias{patternIndex}
\title{Create a Index for a spatial pattern}
\description{Create a Index for a spatial pattern identified from a map
  object, eg \code{link{corField}}, \code{link{mapField}}, or
  \code{link{composite.field}}. The index is computed using a spatial correlation.}
\usage{patternIndex(map,field,anomaly=TRUE)}
\arguments{
    \item{map}{Map describing spatial structure.}
    \item{field}{Field in which the pattern is searched.}
    \item{anomaly}{If true, compute and use the anomalies of the given field.}
 }
 \value{list containing the fields:
   \tabular{rl}{
     index \tab index \cr
     yy    \tab years \cr
     mm    \tab months \cr
     dd    \tab days \cr
   }
   }
\author{R.E. Benestad}}
\examples{
\dontrun{
sst <- retrieve.nc("DNMI_sst.nc")
data(oslo.t2m)
csst <- composite.field(sst,oslo.t2m)
patternIndex(csst,sst,anomaly=FALSE)

lsst<-mapField(sst)
patternIndex(lsst,sst)
}
}
\keyword{manip}