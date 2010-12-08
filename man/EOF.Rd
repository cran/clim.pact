\name{EOF}
\alias{eof}
\alias{EOF}
\alias{Empirical orthogonal Functions}
\alias{PCA}
\alias{ExtEOF}
\alias{lagField}
\alias{EOF2field}
\alias{principal component analysis}
\title{Empirical Orthogonal Functions (EOFs).}
\description{Computes EOFs (a type of principal component
  analysis) for combinations of data sets, typically from
  the NCEP reanalysis and corresponding data from climate models.
  Preprocessing by \code{\link{catFields}} allows for common EOF analysis
  [ref: Benestad (2001), "A comparison between two empirical downscaling
  strategies", \emph{Int. J. Climatology}, \bold{vol 21}, Issue 13,
  pp.1645-1668. DOI 10.1002/joc.703]. and  \code{\link{mixFields}} prepares for
  mixed-field EOF analysis [ref. Bretherton et al. (1992) "An
  Intercomparison of Methods for finding Coupled Patterns in Climate
  Data", \emph{J. Climate}, \bold{vol 5}, 541-560;
  Benestad et al. (2002), "Empirically downscaled temperature scenarios
  for Svalbard", \emph{Atm. Sci. Lett.}, doi.10.1006/asle.2002.0051]. 

  Uncertainty estimates are computed according to North et al. (1982),
  "Sampling Errors in the Estimation of Empirical Orthogonal Functions",
  \emph{Mon. Weather Rev.}, \bold{vol 110}, 699-706. 
  
  NB: This routine may be computer-intensive!
  The computation of the EOFs tends to take some time, especially on
  computers/PCs with little memory (less than 128Mb) and slow processors
  less than 800MHz. 

  See the course notes from Environmental statistics for climate
  researchers \url{http://www.gfi.uib.no/~nilsg/kurs/notes/course.html}
  for a discussion on EOF analysis.

  \code{ExtEOF()} computes extended EOFs.
  
  \code{lagField()} introduces a time lag in the field.

  \code{EOF2field()} reconstructs a field from EOF products. It's a
  useful rountine for filtering out small-scale structure noise if the
  number of EOFs is low - thus including only the few most important
  modes. The function is also useful in conjunction with
  \code{\link{MVR}} to reconstruct the field objects if the regression
  is erformed on EOF objects. The EOF2field function is used when
  reading the gridded data DNMI.sst/DNMI.t2m/DNMI:slp, and is used to
  decompress data stored compressed as EOF-products (note: some detail is lost
  this way). 
  
}

\usage{
EOF(fields,l.wght=TRUE,lc180e=FALSE,direc="data/",
              lon=NULL,lat=NULL,l.stndrd=TRUE,las=1,
              mon=NULL,plot=TRUE,neofs=20,l.rm.ac=TRUE,lsave=TRUE,
              LINPACK=TRUE,silent=FALSE)
ExtEOF(fields,lag=1,mon=NULL,lon=NULL,lat=NULL)
lagField(fields,lag=1)
EOF2field(eof,anomalies=FALSE)
}
\arguments{
    \item{fields}{A field object (eg from \link{retrieve.nc}).}
    \item{l.wght}{'TRUE' applies a geographical weighting.}
    \item{lc180e}{'TRUE' centers the maps on date line (180 deg E).}
    \item{direc}{Directory for the output.}
    \item{lon}{longitudinal region of interest.}
    \item{lat}{latitudinal region of interest.}
    \item{l.stndrd}{Not yet used.}
    \item{las}{Used by \link{filled.contour}, see \link{par}.}
    \item{mon}{Month (1-12) [season (1-4) for daily data] to extract.}
    \item{plot}{'TRUE' plots the results.}
    \item{neofs}{Number of leading EOFs to retain.}
    \item{l.rm.ac}{'TRUE' removes the annual cycle.}
    \item{lsave}{'TRUE' dumps the results to file.}
    \item{LINPACK}{'TRUE': svd; 'FALSE':La.svd}
    \item{silent}{'TRUE': quiet mode}
    \item{lag}{Lag in months or days (using the time unit) introduced to
    the field. Used with \code{\link{mixFields}} to estimate esxtnded
    EOFs.}
    \item{eof}{An 'EOF' object from \code{\link{EOF}}.}
    \item{anomalies}{True if only anomalies are to be shown.}
}
\value{File containing an 'eof.c' object:
  \tabular{ll}{
    EOF \tab EOF patterns. \cr
    W  \tab Eigenvalues. \cr
    dW  \tab uncertainty in Eigenvalues (after North et al., 1982, eq 24: dW approx. W sqrt(2/N)). \cr
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
    clim \tab mean values (climatology) \cr
    attributes \tab Attributes \cr
}
  The data is also saved as files.
}
\author{R.E. Benestad}
\examples{
# Computes a set of mixed-common EOFs (overnight work..). This takes a while...
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
xX.1 <- cat.fields(x.1,X.1,interval.1=c(1958,1998),interval.2=c(1958,2050))
xX.2 <- cat.fields(x.2,X.2,interval.1=c(1958,1998),interval.2=c(1958,2050))
xX <- mix.fields(xX.1,xX.2,mon=1,
                 interval=c(1900,2050))
print("EOF")
eof.c <- EOF(xX.1,mon=1)
eof.mc <- EOF(xX,mon=1)

print("Reconstruct the field from EOFs can filter out small-scale
structure noise:")
EOF2field(xX.1) -> xX.1.filtered

load("data/ceof.Rdata")   # loads the object 'eof'
field <- EOF2field(eof)
}
}
\keyword{spatial}
\keyword{ts}
\keyword{multivariate}

