\name{DNMI.t2m}
\alias{DNMI.t2m}
\alias{DNMI.slp}
\alias{DNMI.sst}
\title{Gridded monthly mean climate data}
\description{Gridded monthly mean 2-meter temperature [T(2m)], sea
level pressure (SLP), and sea surface temperature (SST) analysis for the North
Atlantic region covering the period 1873 -- 1998. The actual data is
filtered though 20 EOFs (see \code{\link{EOF2field}}), and stored as eof products in order to save
space. The lines listed in the example shows how these data have been
generated (NB. old objects with the name \code{eof}, \code{xxx}, \code{climxxx}, and \code{i} 
may be overwritten).}
\usage{data(DNMI.t2m)
       data(DNMI.slp)
       data(DNMI.sst)}
\format{The data is a 'field' object (see \code{\link{retrieve.nc}}).}
\examples{
\dontrun{
library(clim.pact)
t2m <- retrieve.nc("~/data/analysis/DNMI_t2m.nc")
eof <- EOF(t2m)
save(file="clim.pact/data/eof_DNMI_t2m.Rdata",eof)
t2m.eof <- EOF2field(eof)

# Check if the EOF reconstruction reproduces the original field:
newFig()
plotField(t2m,lat=60,lon=10)
plotField(t2m.eof,lat=60,lon=10,add=TRUE,lty=2,col="red")
               # Very similar time series

newFig()
mapField(t2m.eof)
mapField(t2m,add=TRUE,col="red",lty=2,lwd=1)
               # Very similar spatial pattern/contours

sst <- retrieve.nc("~/data/analysis/DNMI_sst.nc")
eof <- EOF(sst)
save(file="clim.pact/data/eof_DNMI_sst.Rdata",eof)

slp <- retrieve.nc("~/data/analysis/DNMI_slp.nc")
eof <- EOF(slp)
save(file="clim.pact/data/eof_DNMI_slp.Rdata",eof)
}
}
\source{The data is also avaialble from \url{http://noserc.met.no}.}
\references{
Benestad, R.E. and Melsom, A. (2002) Is there a link between the
unusually wet autumns autumns in southeastern Norway and SST
anomalies?, Climate Research Vol 23, 67-79.
\\
Melsom E. and Benestad, R.E. :  Is there a link between the unusually
wet autumns in southeastern Norway and SST anomalies?  PDF-version available
on-line from \url{http://www.dams.dk/nmm2002/Melsom\%20and\%20Benestad.pdf}
\\
Benestad, R.E.: Analysis of gridded sea level pressure and 2-meter temperature for 1873-1998 based on UEA and NCEP re-analysis II (RegClim), DNMI KLIMA Report 03/00
  PDF-version available on-line from
  \url{http://noserc.met.no/DS/klima0300.html}. For SST, see
  \url{http://noserc.met.no/DS/dnmisst.html}. The analysis of the SST
  was the same as the T(2m) and SLP.}
\keyword{datasets}
