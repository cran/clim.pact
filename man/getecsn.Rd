\name{getecsn}
\alias{getecsn}
\title{Retrieve daily station record from ECSN set.}
\description{Reads the data from Nordklim available at URL:
  \url{http://eca.knmi.nl/dailydata/predefinedseries.php}.  Also see
  \url{http://www.eumetnet.eu.org/contecsn.html} and \code{\link{getgiss}}.

\bold{THESE DATA CAN BE USED FREELY PROVIDED THAT THE FOLLOWING SOURCE IS ACKNOWLEDGED:}
Klein Tank, A.M.G. and Coauthors (2002) 
'Daily dataset of 20th-century surface air temperature and precipitation series for the European Climate Assessment.'
\emph{Int. J. of Climatol.}, \bold{22}, 1441-1453.
Data and metadata available at \url{http://eca.knmi.nl}

}
\usage{getecsn(location="prompt",param=c("TG","RR"),dataset="blended",
               data.path="data.eca",country=NULL,update = FALSE)
}
\arguments{
  \item{location}{ name of climate station location.}
  \item{param}{ name of element [e.g. avail.elem(), or 't2m', 'rr', 'slp'].}
  \item{dataset}{"blended" (GCOS + EUMETNET), "gcos", or ."eumetnet"}
  \item{data.path}{ name of directory in which the data are stored.}
  \item{country}{Country to search  for stations.}
  \item{update}{For downloading new data.}
}
\value{a <- list of "daily.station.record" class.}
\author{R.E. Benestad}
\examples{
\dontrun{
obs <- getecsn("helsinki")
plotStation(obs)
}
}
\keyword{file}


