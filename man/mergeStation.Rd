\name{mergeStation}
\alias{mergeStation}
\title{Merge climate station series.}
\description{Merges two series from different sources, eg from NACD
  (\code{\link{getnacd}}) and DNMI (\code{\link{getdnmi}}). The code is
  useful for updating long climate series with new observations from a
  different database. The routine cokpares data for overlapping times
  and prints out diagnostics about the two data sets.}
\usage{mergeStation(x.1,x.2,plot=FALSE,print=TRUE,rescale=TRUE)}
\arguments{
  \item{x.1}{1st series.}
  \item{x.2}{2nd series.}
  \item{plot}{'TRUE' plots the overlap.}
  \item{print}{'TRUE' prints diagnostics.}
  \item{rescale}{'TRUE' uses lm to determine amplitude of second sequence otherwise the amplitude of second sequence is not changed.}
  }
\value{A climate station series object.}
\author{R.E. Benestad}}
\examples{
\dontrun{
oslo.1 <- getnacd("OSLO-BLINDERN")
oslo.2 <- getdnmi("oslo")
print(range(oslo.1$yy))
#[1] 1890 1990
print(range(oslo.2$yy))
#[1] 1937 2002
oslo <- mergeStation(oslo.1,oslo.2)
#[1] "Time intervals:"
#[1] 1890 1990
#[1] 1937 2002
#[1] 1937.042 1990.958
#[1] "RMSE:  0.04"
#
#Call:
#lm(formula = y ~ 1 + x, data = ovrlp)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-7.24005 -0.03271  0.01161  0.06006  7.61593
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 0.029044   0.047482   0.612    0.541    
#x           0.993886   0.004866 204.231   <2e-16 ***
#---
#Signif. codes:  0 `***' 0.001 `**' 0.01 `*' 0.05 `.' 0.1 ` ' 1 
#
#Residual standard error: 0.9738 on 644 degrees of freedom
#Multiple R-Squared: 0.9848,	Adjusted R-squared: 0.9848 
#F-statistic: 4.171e+04 on 1 and 644 DF,  p-value: < 2.2e-16

print(range(oslo$yy))
#[1] 1890 2002
}
}
\keyword{manip}
