\name{corEOF}
\alias{field correlation, PCA}
\alias{corEOF}
\title{Field correlation}
\description{Produces maps of spatial correlation patterns from the EOF
  products \deqn{X = {\bf UWV}^T}{X = UWt(V)} according to following formula:
  \deqn{r({\bf X}_{r,t}, \vec{y})=  \frac{
      \sum_k {\bf W}_k {\bf U}_{r,k}{\bf W}_k \sum_t [{\bf V}_{k,t}^T y'_t]}
 {\sqrt{\sum_k{\bf W}_k {\bf U}_{r,k}{\bf W}_k \sum_t({\bf V}_{k,t}^T)^2 
       \times \sum_t(y'_t)^2}}.}{.}

  Reference: Environmental statistics for climate researchers
  \url{http://www.gfi.uib.no/~nilsg/kurs/notes/course.html}

  \bold{NOTE: This routine is not finished yet - still contains some
    bugs, sometimes resulting in absolute correlation values greater
    than unity.} 
}
\usage{corEOF(x,y,lsig.mask=TRUE,sig.lev=0.05,neofs=20,
                 lty=1,col="black",lwd=1)}
\arguments{
    \item{x}{An EOF.}
    \item{y}{A station series.}
    \item{lsig.mask}{Mask out values that are not statistically significant.}
    \item{sig.lev}{Level of significanse.}
    \item{neofs}{Number of modes to include.}
    \item{lty}{Contour line type.}
    \item{col}{Contour line colour.}
    \item{lwd}{Contour line width.}
 }
\value{A map object.}
\author{R.E. Benestad}}
\examples{
data(oslo.t2m)
data(eof.slp)
corEOF(eof.slp,oslo.t2m)
}
\keyword{manip}