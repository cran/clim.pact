\name{dist2norm}
\alias{dist2norm}
\alias{norm2dist}
\title{Transform a series to a normally distributed series.}
\description{Used for DS on non-Gaussian data. Tailored for the use with
  daily precipitation. \code{dist2norm} transforms (maps) data to normal
distribution and \code{norm2dist} carries out the reverse operation. The
mapping/transformation is done through one-to-one matching of the
empirical distribution functions (e.d.f.) of the data and a Gaussian
distribution.}
\usage{
  dist2norm(x,plot=FALSE,exclude=NULL,sd=1,mean=0,force.zero=TRUE)
  norm2dist(x,plot=FALSE,exclude=NULL,sd=1,mean=0,force.zero=TRUE)
}
\arguments{
    \item{x}{vector or a list-object returned by \code{dist2norm}.}
    \item{plot}{T -> plot for illustrating the mapping.}
    \item{exclude}{Values (e.g. 0) to exclude for the transform.}
    \item{sd}{The standard deviation of the transformed quantities.}
    \item{mean}{The mean value of the transformed quantities.}
    \item{force.zero}{Force the e.d.f. to start from 0.}
 }
\value{a list object of class 'dist2norm' containing a vector and the empirical transform functions.}
\author{R.E. Benestad}
\examples{
data(oslo.dm)
x <- dist2norm(oslo.dm$precip)
plot(x$xT)

y <- norm2dist(x,plot=TRUE)
}
\keyword{manip}


