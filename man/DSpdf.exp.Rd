\name{Downscale PDF}
\alias{DSpdf.exp}
\alias{exp.law1}
\alias{exp.par}
\alias{CDFtransfer}
\alias{empiricalRanking}
\title{Downscale eponential PDF.}
\description{
  Functions handy for downscaling PFDs.

  DSpdf.exp() is used to downscale expeoential distributions for daily
  rainfall according to Benestad, R.E., C. Achberger, &  E. Fernandez
  (2005) Empirical-statistical downscaling of distribution functions for
  daily precipitation, met-no report, 12-2005, Climate, pp.43
  (\url{http://met.no/english/r_and_d_activities/publications/2005/12_2005/abstract_12_2005.html}). The dataset exp.law1 is used for this downscaling. 

  CDFtransfer() is a function that facilitates local quantile
  transformation.

  empiricalRanking() uses a formula for estimating the cumulative probability P corresponding to rank m
 References: Jenkinson, A.F., 1977, U.K. Met.Office Synoptic
 Clim. Branch Memo 58;  Beard, L.R., 1943, Trans. Amer. Meteor. Soc. Civ. Eng., 108, 1110-1160; Chegodaev, N.N., 1953 (in Russian) State Rail Transport Publishing House; Folland, C. and Anderson, C. (2002), J. Clim. 15, 2954-2960, equation (1)
}
\usage{
DSpdf.exp(obs=NULL,dT=0,dP=0,plot=TRUE,year=NULL,month=NULL)
data(exp.law1)
CDFtransfer(Y,CDF.2,CDF.1=NULL,method="empiricalRanking",plot=FALSE,silent=FALSE,smooth=TRUE)
empiricalRanking(x)
}
\arguments{
  \item{x}{A field object.}
  \item{plot}{TRUE: plot.}
  \item{obs}{daily station object}  
  \item{dT}{Projected change in annual mean T(2m), e.g. from objDS (unit
    deg C).}  
  \item{dP}{Projected change in annual mean precipitation, e.g. from
    objDS (units: mm/day!).}
  \item{year}{Year to extract}
  \item{month}{Month to extract}
  \item{Y}{A data series or daily station object.}
  \item{CDF.1}{If provided, assume this CDF to describe the cumulative
    distribution of Y.}
  \item{CDF.2}{cumulative distribution predicted for the future.}
  \item{method}{Method for estimating emprical distribution function (EDF).}
  \item{silent}{TRUE: less clutter on the screen.}
  \item{smooth}{Uses a 9th-order polynomial fit to provide a smoother
    fit for the local quantile transformation.}
 }
\value{A field object.}
\author{R.E. Benestad}
\examples{
data(exp.law1)
data(oslo.dm)
a<-DSpdf.exp(oslo.dm,dT=3,dP=1)
F1<- list(x=a$x,P=a$Fx.obs)
F2<- list(x=a$x,P=a$Fx.chg) 
y<-CDFtransfer(Y=oslo.dm$precip,CDF.2=F2,CDF.1=F1,plot=TRUE)
plot(oslo.dm$precip,y)
}
\keyword{hplot}
