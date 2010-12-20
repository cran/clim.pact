\name{mixFields}
\alias{mixFields}
\alias{combineFields}
\alias{combineEOFs}
\title{mixFields}
\description{Mix fields by combining two different gridded sets of
  observation. Observations/data for representing values at n different
  locations at a given time (t) can be described in terms of a
  vector
  \deqn{\vec{x}(t) = [x_1, x_2, ... x_n].}{x(t)=[x_1, x_2,
    ... x_n].}
  Two different sets of observations can be represented by two vectors
  y and z of lengths n and m respectively. In mix.fields, the
  information in these two data sets are combined combining the two vectors:
  \deqn{\vec{x}(t) = [\vec{y}(t),\vec{z}(t)] =
    [y_1, y_2, ... y_n, z_1, z_2, ... z_m].}{
    x(t)=[y(t), z(t)]= [y_1, y_2, ... y_n, z_1, z_2, ... z_m].}
  The length of the final vector of the mixed field is the sum of the
  lengths of the two respective vectors. The two data sets do not have
  to be on the same grid.

  reference: Bretherton et al. (1992) "An
  Intercomparison of Methods for finding Coupled Patterns in Climate
  Data", \emph{J. Climate}, \bold{vol 5}, 541-560.

  The output from \code{mixFields} can be used in \code{\link{EOF}}
  to compute mixed-common EOFs which subsecuently can be used as
  predictors in \code{\link{DS}} in order to downscale climate scenarios
  (Benestad et al. (2002), "Empirically downscaled temperature scenarios
  for Svalbard", \emph{Atm. Sci. Lett.}, doi:10.1006/asle.2002.0050).

  The functions \code{combineFields} and \code{combineEOFs} bundle
  together fields of different parameters, and allows downscaling of
  multiple variables, but without performing a mixed-EOF analysis. Thus,
  allows ESD with a multiple regression analysis involving both e.g. SLP
  and T(2m) as input variables as in traditional multiple regression:

  \deqn{\hat{y}(t) = alpha_0 + alpha_1 x(t) + alptha_2 + z(t) +
    ....}{y(t)=a0 + a1 x(t) + a2 + z(t) + ...}  
 
}  
\usage{
mixFields(field.1,field.2,mon=NULL,interval=NULL)
combineFields(field.1,field.2)
combineEOFs(eof.1,eof.2)
}
\arguments{
    \item{field.1}{A 'field.object'.}
    \item{field.2}{A 'field.object'.}
    \item{eof.1}{An 'EOF.object'.}
    \item{eof.2}{An 'EOF.object'.}
    \item{mon}{Calendar month to extract.}
    \item{interval}{Time interval to extract.}
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


}
}
\keyword{models}
