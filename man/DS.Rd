\name{DS}
\alias{climate analysis}
\alias{DS}
\alias{downscaling}
\title{Downscaling of monthly means}
\description{Identifies statistical relationships between large-scale
  spatial climate patterns and local climate variations for monthly and
  daily data series. Calibrates a linear regression model using
  step-wise screening and common EOFs (\code{\link{EOF}}) as basis
  functions. Evaluates the statistical relationship. Predicts local
  climate parameter from predictor fields. Works with ordinary EOFs,
  common EOFs (\code{\link{catFields}}) and mixed-common EOFs
  (\code{\link{mixFields}}).
  The rationale for using mixed-common EOFs is that the coupled
  structures described by the mixed-field EOFs may have a more physical
  meaning than EOFs of single fields [Benestad et al. (2002),
  "Empirically downscaled temperature scenarios for Svalbard",
  \emph{Atm. Sci. Lett.}, doi.10.1006/asle.2002.0051].

  The downscaling analysis returns a time series representing the local
  climate, patterns of large-scale anomalies associated with this,
  ANOVA, and analysis of residuals. Care must be taken when using this
  routine to infer local scenarios: check the R2 and p-values to check
  wether the calibration yielded an appropriate model. It is also
  important to examine the spatial structures of the large-scale
  anomalies assocaiated with the variations in the local climate: do
  these patterns make physical sense? Experiment with both single and
  mixed fields. It is also a good idea to check
  whether there are any structure in the residuals: if so, then a linear
  model for the relationship between the large and small-scale
  structures may not be appropriate. It is furthermore important to
  experiment with predictors covering different regions [ref: Benestad
  (2001), "A comparison between two empirical downscaling strategies",
  \emph{Int. J. Climatology}, \bold{vol 21}, Issue 13, pp.1645-1668. DOI
  10.1002/joc.703]. 

  The function \code{ds()} is a generic routine which in principle works for
  when there is any real statistical relationship between the predictor
  and predictand. The predictand is therefore not limited to a climate
  variable, but may also be any quantity affected by the regional
  climate. \emph{It is important to stress that the downscaling model
  must reflect a well-understood (physical) relationship.} 

  The trend-estimation uses regression to fit a 5th-order polynomial (in
  time) to fit the observed time series. The rate-of-change is estimated
  by taking the time-derivative of this equation. If \deqn{y= c_0 + c_1
  x + c_2 x^2 + c_3 x^3 + c_4 x^4 + c_5 x^5,}{y= c0 + c1 x + c2 x^2 + c3
  x^3 + c4 x^4 + c5 x^5,} where x is the time, then the rate-of-change
is: \deqn{y= c_1 + 2 c_2 x + 3 c_3 x^2 +  4 c_4 x^3 + 5 c_5 x^4.}{y= c1
  + 2 c2 x + 3 c3 x^2 + 4 c4 x^3 + 5 c5 x^5.} [ref: Benestad
  (2002), What can present climate models tell us about climate
  change?, \emph{Climatic Change}, accepted.] 
  
  The routine uses a step-wise regression (step) using the leading
  EOFs. The calibration is by default carried out on de-trended data
  [ref: Benestad (2001), "The cause of warming over Norway in the
  ECHAM4/OPYC3 GHG integration", \emph{Int. J. Clim.}, 15 March,
  \bold{vol 21}, p.371-387.].

  The downscaled scenario is saved in a text file in the output
  directory (default: 'output').

  The course notes from Environmental statistics for climate researchers
  \url{http://www.gfi.uib.no/~nilsg/kurs/notes/course.html} is a useful
  reference to statistical modelling and regression. 
}
\usage{DS(dat,preds,mon=NULL,direc="output/",cal.id=NULL,
              ldetrnd=TRUE,i.eofs=seq(1,8,by=1),ex.tag="",
              method="lm",plot=TRUE,leps=FALSE,param="t2m",
              plot.res=FALSE,plot.rate=FALSE,xtr.args="")}

\arguments{
    \item{dat}{ A climate.station object (\code{\link{station.obj}} or 
        \code{\link{station.obj.dm}}). 
        [e.g. from \code{getnacd}, \code{getnordklim} or \code{station.obj}].}
    \item{preds}{ The predictor \code{\link{EOF}}.}
    \item{mon}{ month or season to downscale, this is automatically 
        changes if predictor only contains a different
	month (this is normally a redundant feature).}
    \item{direc}{ name of directory inwhich the output is dumped
          (e.g. figures, tables).}
    \item{cal.id}{ID tag used for calibration. By default use the 
	first field (\code{\link{catFields}}) for calibration.}
    \item{ldetrnd}{ F for no detrending; T for removing linear trends
	before model calibration.}
    \item{i.eofs}{ select which EOFs to include in the setp-wise screening.}
    \item{ex.tag}{Extra labelling tag for file names for experiments.}
    \item{method}{Sets the method to use for regression.}
    \item{plot}{'TRUE' produces figures.}
    \item{leps}{'TRUE' produces EPS figures (files).}
    \item{param}{Name of parameter (for plot labels).}
    \item{plot.res}{'TRUE' shows statistics for residuals.}
    \item{plot.rate}{'TRUE' shows analysis of rate-of-change.}
    \item{xtr.args}{Extra/additional arguments in the formula.}
}

\value{A 'ds' object -  a list of following elements:
  \tabular{rl}{
    X.1 .. X.n \tab 1..nth predictor pattern for n fields
    (\link{mixFields}). \cr
    lon.1 .. \tab Longitude coordinate of spatial fields (a vector). \cr
    lat.1 .. \tab Latitude coordinate of spatial fields (a vector). \cr
    n.fld \tab Number of fields (different types of predictors,
    \link{mixFields}). \cr
    unit \tab Unit of quantity in station series. \cr
    pred.name \tab Name of predictor. \cr
    lon.loc \tab Longitude of predictand location. \cr
    lat.loc \tab Latitude of predictand location  \cr
    yy.gcm \tab Years corresponding to scenario (GCM). \cr
    mm.gcm \tab Months corresponding to scenario (GCM). \cr
    dd.gcm \tab Days corresponding to scenario (GCM). \cr
    yy.cal \tab Years corresponding to observation (Calibration). \cr
    mm.cal \tab Months corresponding to observation (Calibration). \cr
    dd.cal \tab Days corresponding to observation (Calibration). \cr
    yy.o \tab Years corresponding to station series (obs.). \cr
    mm.o \tab Months corresponding to station series (obs.). \cr
    dd.o \tab Days corresponding to station series (obs.). \cr
    rate.ds \tab Estimated linear rate of change of downscaled scenario. \cr
    rate.err \tab Error estimate for rate.ds. \cr
    gcm.trnd.p \tab P-value of linear trend in downscaled scenario.\cr
    fit.p \tab ANOVA p-value for fit between large-scale and small-scale \cr
          \tab variability(from regression analysis). \cr
    fit.r2 \tab ANOVA R2 for fit between large-scale and small-scale \cr
           \tab variability (from regression analysis). \cr
    pre.gcm \tab The downscaled scenario (a vector). \cr
    pre.y \tab The downscaled results using the calibration data. \cr
    location \tab Nsme of location of predictor. \cr
    gcm.stat \tab ANOVA of linear trend fit to scrnario. \cr
    month \tab Month of study (0-> all months).\cr
    v.name \tab Name of downscaled element. \cr
    region \tab Region used for downscaling. \cr 
    pre.fit \tab Linear fit to prediction (downscaled scenario) (a vector). \cr
    pre.p.fit \tab Polynomial fit to the downscaled scenario.\cr  
    tr.est.p.fit \tab Rate of change derived from a fifth-order
    polynomial  \cr
    \tab trend-fit to prediction (downscaled scenario) (a vector).\cr
    id.1, id.2 \tab IDs labelling which data was used for calibration (id.1).
  }
}
\author{R.E. Benestad}}
\examples{
library(clim.pact)
data("oslo.t2m")
data("eof.mc")
a<-DS(dat=oslo.t2m,preds=eof.mc,plot=FALSE)
}
\keyword{models}
\keyword{multivariate}
\keyword{ts}
\keyword{spatial}