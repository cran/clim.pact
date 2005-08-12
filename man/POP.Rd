\name{POP}
\alias{Principal Oscillating Patterns}
\alias{POP}
\alias{plotPOP}
\title{Principal Oscillating Patterns}
\description{After von Storch & Zwiers (1999), Statistical Analysis in Climate Research, p. 338.
%\deqn{X_{t+1} = A X_t + \eta}{X(t+1) = A X(t) + noise}.
}
\usage{
POP(x,plot=TRUE,main="POP analysis",sub="",test=FALSE,i.eofs=1:8,LINPACK=TRUE, mode = 1)
plotPOP(pop,mode = 1, main = "POP analysis", sub = "")
}
\arguments{
    \item{x}{A field or an eof object.}
    \item{plot}{Flag: plot the diagnostics.}
    \item{test}{Flag: test by reconstructing one series (leading EOF or a grid-box series).}
    \item{i.eofs}{Which EOFs to include (only when the input is given as eof objects).}
    \item{LINPACK}{'TRUE': svd; 'FALSE':La.svd}
    \item{main}{main title (see \code{link{plot}}).}
    \item{sub}{subtitle (see \code{link{plot}}).}
    \item{mode}{POP mode to display}
    \item{pop}{a POP object.}
 }
\value{A POP object: an  \code{link{eigen}} object with additional variables. maps (map of the absolute of POP patterns), pop.Im (the imaginary POP pattern), pop.Re(the real POP pattern), lon (longitude coordinates), and lat (latitude coordinates).}
\author{R.E. Benestad}}
\examples{
\dontrun{
data(DNMI.t2m)
eof.1 <- EOF(DNMI.t2m,mon=1)
pop <- POP(eof.1)
}
}
\keyword{manip}
