\name{EOF2field}
\alias{EOF2field}
\title{Reconstructs a field from EOF products}
\description{Reconstructs a field from \code{\link{EOF} products}.
}
\usage{EOF2field(eof,anomalies=FALSE)}

\arguments{
    \item{eof}{An 'EOF' object from \code{\link{EOF}}.}
    \item{anomalies}{True if only anomalies are to be shown.}
}

\value{A 'field' object (see \code{\link{retrieve.nc}}).}
}
\author{R.E. Benestad}
\examples{
\dontrun{
load("data/ceof.Rdata")   # loads the object 'eof'
field <- EOF2field(eof)
}
}
\keyword{models}
\keyword{multivariate}
\keyword{ts}
\keyword{spatial}
