\name{optint}
\alias{optimum interpolation}
\alias{optint}
\title{Optimum interpolation}
\description{
  Solves for w in:
  \deqn{\sum_j M_{ij} w_ik = \langle pi_j pi_k \rangle }{sum_j(M_ij
    w_ik) = <pi_j pi_k>. }

  Then computes a correction: \deqn{r = w_i q}{r = w_i q,}
  where q = gridded observations - first guess.

  Reference: Reynolds and Smith (1994), Improved global sea surface
  temperature analysis using optimum interpolation, \emph{J. Climate}, 
  \bold{vol 7}, 929-948.

  \bold{NB: The testing of this routine is not yet finished - there may
    still be some bugs in this function.}
   }
\usage{optint(lon,lat,obs,lon.grd,lat.grd,fguess,eps,
              lambda=50,M=NULL,piipij=NULL,w=NULL,tim=NULL,date=NULL)}
\arguments{
    \item{lon}{A vector of longitude coordinates.}
    \item{lat}{A vector of latitude coordinates.}
    \item{obs}{A vector of observations.}
    \item{lon.grd}{A vector of longitude coordinates for gridded values.}
    \item{lat.grd}{A vector of latitude coordinates for gridded values.}
    \item{fguess}{First guess values on given grid (\code{lon.grd,lat.grd}).}
    \item{eps}{Error estimates on given grid (\code{lon.grd,lat.grd}).}
    \item{lambda}{Correlation length.}
    \item{M}{Interpolation model matrix.}
    \item{piipij}{Interpolation weights.}
    \item{w}{Interpolation weights..}
    \item{tim}{Time index.}
    \item{date}{Date stamp.}
 }
\value{A map object with additional elements: \code{M, piipij, w}}
\author{R.E. Benestad}}
\examples{
# This example takes a little while to run..
\dontrun{
lon <- c(12.02, 10.47,  9.28, 10.78, 10.71,  9.12,  8.20,  9.52,  8.18,  8.52,
          8.07,  6.15,  5.63,  6.37,  5.92,  5.33,  6.50,  7.42,  6.17, 10.00,
         10.02, 11.12, 10.93, 14.02, 13.98, 14.43, 15.42, 17.80, 18.53, 18.93,
         21.02, 21.88, 23.35, 23.52, 31.08, 15.47, 11.93)
lat <- c(60.61, 61.10, 62.12, 59.38, 59.95, 60.58, 60.52, 59.90, 59.83, 59.03,
         58.20, 58.88, 58.88, 59.65, 60.40, 60.38, 60.65, 62.23, 62.57,
         63.20, 62.82, 63.20, 63.47, 65.52, 66.82, 67.27, 68.70, 68.73,
         69.06, 69.65, 69.73, 69.83, 69.98, 69.33, 70.37, 78.25, 78.92)
val <- c(1.74, 1.43, 0.98, 1.31, 1.10, 1.39, 1.07, 1.33,-0.04, 1.07,
         1.07, 0.22, 0.89, 0.98, 0.77, 0.81, 1.43, 0.30, 0.69,-0.17,
         1.17, 0.32, 1.55, 1.65, 0.75, 0.78, 0.56, 0.04, 0.57, 0.52,
         0.51, 0.51, 0.96, 0.32, 0.35,-0.36, 0.00) 
library(akima)
lon.grd <- seq(min(lon),max(lon),length=20)
lat.grd <- seq(min(lat),max(lat),length=20)
fguess <- interp(lon,lat,val,lon.grd,lat.grd)$z
fguess[is.na(fguess)]<-mean(fguess,na.rm=TRUE)
eps <- abs(fguess)*0 + 0.1
opt.int <- optint(lon,lat,val,lon.grd,lat.grd,fguess,eps)
image(opt.int$lon,opt.int$lat,opt.int$map)
contour(opt.int$lon,opt.int$lat,opt.int$map,lwd=1,add=TRUE)
contour(lon.grd,lat.grd,fguess,lty=2,lwd=2,col="darkblue",add=TRUE)
addland()
}
}
\keyword{manip}
