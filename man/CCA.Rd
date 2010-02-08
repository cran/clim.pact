\name{CCA}
\alias{Canonical correlation analysis}
\alias{CCA}
\alias{testCCA}
\alias{Psi}
\alias{predictCCA}
\alias{plotCCA}
\alias{stations2field}
\alias{check.repeat}
\title{Canonical correlation analysis}
\description{
Applies a canonical correlation analysis (CCA) to two data sets. The CCA here
can be carried out in two different ways: i) an \code{\link{svd}} based approach (after 
Bretherton et al. (1992), J. Clim. Vol 5, p. 541, also documented in Benestad (1998): 
"Evaluation of Seasonal Forecast Potential for Norwegian  Land Temperatures and Precipitation using CCA", 
DNMI KLIMA Report 23/98 at \url{http://met.no/english/r_and_d_activities/publications/1998.html}) or ii) a 
covariance-eigenvalue approach (after Wilks, 1995, "Statistical methods in the Atmospheric 
Sciences", Academic Press, p. 401). 

The analysis can also be applied to either EOFs or fields.

Note: the analysis has sometimes been somewhat unstable, returning inconsistent results. The recommendation is to use EOFs and SVD option.

The CCA analysis can be used to develope statistical models according to:

\deqn{Y = \Psi X}{Y = Psi X}

Where Y is the predictand and X the predictor.
\code{plotCCA} plots the CCA results, \code{testCCA} is for code
verification, and \code{Psi} returns the matrix \deqn{\Psi}{Psi}.

\code{stations2field} turns a group of station objects into a field by
the means of a simple and crude interpolation/gridding. \code{check.repeat} is a quality-control
function that eliminates repeated years in the station objects.
}
\usage{
CCA(x1,x2,SVD=TRUE,plot=TRUE,main="CCA",sub="",
    test=FALSE,i.eofs=1:8,LINPACK=TRUE)
plotCCA(cca)
testCCA(method="CCA",reconstr=FALSE,mode=1,test=TRUE,
        LINPACK=TRUE,SVD=TRUE,n.pc=4,synthetic=TRUE)
Psi(cca)
predictCCA(Psi,X)
stations2field(data.set=c("narp"),ele=101,obj.type="monthly.field.object",
               plot=TRUE,silent=FALSE,intrp.method="interpp",
               interpolation.option="simple")
check.repeat(x)

}
\arguments{
    \item{cca}{A CCA object.}
    \item{x1}{A field or an eof object.}
    \item{x2}{A field or an eof object.}
    \item{SVD}{Flag: determine which approach to use: SVD or eigenfunction-based alogithm.}
    \item{plot}{Flag: plot the diagnostics.}
    \item{test}{Flag: test by reconstructing one series (leading EOF or a grid-box series).}
    \item{i.eofs}{Which EOFs to include (only when the input is given as eof objects).}
    \item{LINPACK}{'TRUE': svd; 'FALSE':La.svd}
    \item{main}{main title (see \code{link{plot}}).}
    \item{sub}{subtitle (see \code{link{plot}}).}
    \item{method}{Which method to test: CCA or MVR.}
    \item{reconstr}{For the test-reconstruction of fields.}
    \item{mode}{Test for a particular EOF pattern/mode - the other modes are randomized.}
    \item{n.pc}{Number of principal components to include.}
    \item{synthetic}{Construct artificial test data from a random number
      generator and cosine series.}
    \item{Psi}{CCA-based prediction model (a matrix): Y = Psi X}
    \item{X}{A field object used as predictor: : Y = Psi X}
    \item{data.set}{Strings (eg, "narp","nordklim",or "nacd") if data
      has been installed, or list of station objects (and nothing else)}.
    \item{ele}{element code - see \code{\link{getnordklim}}.}
    \item{obj.type}{Object type, e.g "monthly.field.object"}
    \item{silent}{FALG - verbose or not.}
    \item{intrp.method}{Method fro gridding the data.}
    \item{interpolation.option}{Either "simple","distance", or
      "test". In current version, it's only "simple" which seems to work.}
    \item{x}{station object.}
 }
\value{A CCA object: a list containing a.m, b.m, u.k, v.k, and r, describing the Canonical Correlation variates, patterns and correlations. a.m and b.m are the patterns and u.k and v.k the vectors (time evolution).}
\author{R.E. Benestad}
\examples{
\dontrun{
data(DNMI.t2m)
data(DNMI.slp)
eof.1 <- EOF(DNMI.t2m,mon=1)
eof.2 <- EOF(DNMI.slp,mon=1)
cca <- CCA(eof.1,eof.2)
# Testing routine:
testCCA()


}
}
\keyword{manip}
