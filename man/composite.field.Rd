\name{composite.field}
\alias{composite.field}
\alias{compositeField}
\title{Composite maps}
\description{Produce composites of maps based on station observations,
  weights, or years. For a station object, assign +/- 1 weights
  according to whether the values stored in \code{val} (monthly station
  objects)  are greater than mean + sd or less than
  mean - sd.}
\usage{
composite.field(x,y,lsig.mask=TRUE,sig.lev=0.05,s=0.42,mon=NULL,
                      lty=1,col="black",lwd=1)
compositeField(x,y,lsig.mask=TRUE,sig.lev=0.05,s=0.42,mon=NULL,
                      lty=1,col="black",lwd=1)
}
\arguments{
    \item{x}{A field object.}
    \item{y}{A station object, a vector consisting of [-1,0,+1] or a
      vector consisting of the years (negative values are used for
      negative phase, eg. c(1991,1993,1998,-1961,-1963,-1994).}
    \item{lsig.mask}{FLAG: mask the regions not statistically significant.}
    \item{sig.lev}{Level of significance.}
    \item{s}{Threshold for defining whether \code{y} is high or low = s *
      sd(y).}
    \item{mon}{Month to analyse.}
    \item{lty}{Contour line type.}
    \item{col}{Contour line colour.}
    \item{lwd}{Contour line width.}
 }
\value{A map object.}
\author{R.E. Benestad}}
\examples{
\dontrun{
slp <- retrieve.nc("ncep_slp.nc")
data(oslo.t2m)
composite.field(slp,oslo.t2m)
}
}
\keyword{manip}
