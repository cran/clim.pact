\name{mergeEOF}
\alias{mergeEOF}
\title{Merge EOFs}
\description{Merges two EOF objects by making the first one resembling
the second (a bit similar to \code{\link{mergeStation}}). Can be used
for reconstructing fields.} 
\usage{mergeEOF((eof1,eof2,plot=TRUE,silent=FALSE,method="lm",
               match="time",cut.off=8,adjust=TRUE))}
\arguments{
    \item{eof1}{The EOF object (see \code{\link{EOF}}) for the earliest data (assuming this
    data set contains most errors).}
    \item{eof2}{The EOF object (see \code{\link{EOF}}) for the most recent data (assuming this
    data set contains less errors)}
    \item{plot}{if 'TRUE' then plot diagnostics}
    \item{silent}{if 'TRUE' then do not print out diagnostics.}
    \item{method}{"lm" or "matrix projection" (G. Strang (1988),
    "Linear algrebra and its applications", Hartcourt Brace & Company,
    3rd ed.(p.147))}
    \item{match}{"time" or "space"}
    \item{cut.off}{How many EOFs to include}
    \item{adjust}{TRUE: apply \code{\link{adjust.eof}} to results.}
  }
\value{An 'eof' object.}
\author{R.E. Benestad}}
\examples{
\dontrun{
data(DNMI.slp)
NCEP.slp<-retrieve.nc("~/data/ncep/slp.mon.mean.nc",
                    x.rng=c(-90,50),y.rng=c(0,75))

#--------------------------------------------------
# Need to fix some details of the NCEP.slp object (only for 
# proper 'housekeeping')

NCEP.slp$dd[] <- 15
attr(NCEP.slp$tim,unit') <- "month"
class(NCEP.slp) <- c("field","monthly.field.object")
#---------------------------------------------------

eof1 <- EOF(DNMI.slp,mon=1)
eof2 <- EOF(NCEP.slp,mon=1)
eof <- mergeEOF(eof1,eof2)
}
}
\keyword{manip}
