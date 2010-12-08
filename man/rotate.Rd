\name{rotate}
\alias{rotate}
\title{Rotate spherical coordinates}
\description{The function uses the 3-dimensional space to find vectors
  describing the locations on a sphere that correspond to a given
  longitude and latitude. Two reference points are given, one for the
  new 'north pole': 
  \eqn{\vec{r}_0 = [\cos(\phi_0)\cos(\theta_0),\sin(\phi_0),\cos(\phi_0)\sin(\theta_0)]}{r0 = [cos(phi0)cos(theta0),sin(phi0),cos(phi0)\sin(theta0)]} and
  one for a point on the same meridian \eqn{\vec{r}_x}{rx}, but closer to the
  equator. The vector corresponing to each of the longitude-latitude point is given by:
  \eqn{\vec{r}_i = [\cos(\phi_i)\cos(\theta_i),\sin(\phi_i),\cos(\phi_i)\sin(\theta_i)]}{r = [cos(phi)cos(theta),sin(phi),cos(phi)sin(theta)]}.

  The default method is based on "Cayley-Klein"
  (\url{http://en.wikipedia.org/wiki/Rotation_representation_(mathematics)}),
  but an optional approach uses inner products of 3D-vectors where the new spherical coordinates are computed as follows: new latitude
  (\code{phi}) is the angle between the the vectors \eqn{\vec{r}_0}{r0}
    and \eqn{\vec{r}_i}{r[i]} and estimated from the inner products (\url{http://en.wikipedia.org/wiki/Inner_product_space}). The new longitude (\code{theta}) is the
      angle between the the vectors \eqn{\vec{r}_i-\vec{r}_0}{r[i]-r0} 
    and \eqn{\vec{r}_x-\vec{r}_0}{rx-r0}. The angles are estimated by
    taking the inner-product and the arc-cosine.
}
\usage{rotate(lons,lats,lon.0=NULL,lat.0=NULL, method = "Cayley-Klein", test = TRUE)}
\arguments{
    \item{lons}{Longitude coordinates to be rotated.}
    \item{lats}{Latitude coordinates to be rotated.}
    \item{lon.0}{Longitude coordinate of point corresponding to new 'north'pole.}
    \item{lat.0}{Latitude coordinate of point corresponding to new 'north'pole.}
    \item{method}{Method to use for rotation, e.g. "Cayley-Klein" or "old"}
    \item{test}{Flag for testing the function.}
 }
\value{list containing \code{phi} and \code{theta}, the new spherical coordinates.}
\author{R.E. Benestad}
\examples{
data(addland1)
rot <- rotate(lon.cont,lat.cont,lon.0=0,lat.0=0)
}
\keyword{manip}
