\name{MVR}
\alias{Multivariate regression analysis}
\alias{MVR}
\title{Multivariate regression analysis}
\description{Applies a multivariate regression (MVR) analysis to two data sets. The MVR here
is based on the projection to obtain a least squares approximation and uses the formula of 
Strang (1988) "Linear Algebra and its applications", Harcourt Brace and Company, p. 156. 
The method is also documented in Benestad (1999) "MVR applied to Statistical Downscaling for 
prediction of Monthly Mean Land Surface Temperatures: Model Documentation", DNMI KLIMA 
Report 02/99 at \url{http://met.no/english/r_and_d_activities/publications/1999.html}.

For the expression \deqn{Ax = b}{Ax=b}, then the projection of b onto the columns space of A through
\deqn{p= A(A^T A)^{-1}B^T b}{p= A * inv[t(A) * A] * t(B) * b}.

Here \code{x} is called the predictor, \code{b} is the prediction of the
equation, and the data - called the predictand - used for calibration is
\code{y}.

The predictor and predictand arguments can be both EOF-objects or field
objects. The former requires much less memory.

Also see \code{\link{CCA}}, \code{\link{EOF}} and \code{\link{EOF2field}}.
}

\usage{
MVR(x,y,plot=TRUE,main="Multivariate regression",sub="",test=FALSE,
    i.eofs=1:8,LINPACK=TRUE, SVD=TRUE)
}
\arguments{
    \item{x}{A field or an eof object. The predictor}
    \item{y}{A field or an eof object. The predictand}
    \item{plot}{Flag: plot the diagnostics.}
    \item{test}{Flag: test by reconstructing one series (leading EOF or a grid-box series).}
    \item{i.eofs}{Which EOFs to include (only when the input is given as eof objects).}
    \item{LINPACK}{'TRUE': svd; 'FALSE':La.svd}
    \item{main}{main title (see \code{link{plot}}).}
    \item{sub}{subtitle (see \code{link{plot}}).}
    \item{SVD}{Flag: determine which approach to use: SVD or eigenfunction-based alogithm.}
 }
\value{A MVR object that is similar to a field or EOF object (inherits the object type, with an additional "MVR" label) with the projection (dat), but with additional fields such as the weights (psi, which is a map object) [dat represents p and psi represents x.hat in Strang (1988)].}
\author{R.E. Benestad}
\examples{
\dontrun{
data(DNMI.t2m)
data(DNMI.slp)

# MVR with fields
DNMI.t2m.jan <- catFields(DNMI.t2m,mon=1)
DNMI.slp.jan <- catFields(DNMI.slp,mon=1)
mvr.field <- MVR(DNMI.t2m.jan,DNMI.slp.jan)

# MVR with EOFs
eof.1 <- EOF(DNMI.t2m,mon=1)
eof.2 <- EOF(DNMI.slp,mon=1)
mvr.eof <- MVR(eof.1,eof.2)

DNMI.t2m.jan <- catFields(DNMI.t2m,mon=1)
DNMI.slp.jan <- catFields(DNMI.slp,mon=1)
mvr.field.1 <- MVR(DNMI.t2m.jan,DNMI.slp.jan,LINPACK=FALSE)
mvr.field.2 <- catFields(MVR(DNMI.t2m,DNMI.slp),mon=1)
mvr.field.3 <- MVR(DNMI.t2m.jan,DNMI.slp.jan)
EOF(mvr.field.3)-> eof3
EOF(mvr.field.2)-> eof2
EOF(mvr.field.1)-> eof1
plotEOF(eof1)
plotEOF(eof2)
plotEOF(eof3)
}
}
\keyword{manip}
