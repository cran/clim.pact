\name{corField}
\alias{field correlation}
\alias{corField}
\title{Field correlation}
\description{Produces maps of spatial correlation patterns.}
\usage{corField(x,y,lsig.mask=TRUE,sig.lev=0.05,mon=NULL,
                 lty=1,col="black",lwd=1,main=NULL,z.levs=NULL,my.col=NULL)}
\arguments{
    \item{x}{A field object.}
    \item{y}{A station series or a field object.}
    \item{lsig.mask}{Mask out values that are not statistically significant.}
    \item{sig.lev}{Level of significanse.}
    \item{mon}{Month to analyse.}
    \item{lty}{Contour line type.}
    \item{col}{Contour line colour.}
    \item{lwd}{Contour line width.}
    \item{main}{Title of plot, same as 'main' in plot().}
    \item{z.levs}{contour levels.}
    \item{my.col}{colour template (see \code{link{rgb}}).}
 }
\value{A map object.}
\author{R.E. Benestad}}
\examples{
\dontrun{
slp <- retrieve.nc("ncep_slp.nc")
data(oslo.t2m)
corField(slp,oslo.t2m)
}
}
\keyword{manip}
