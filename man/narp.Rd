\name{narp}
\alias{narp}
\alias{narp.meta}
\alias{nacd.meta}
\alias{nordklim.meta}
\title{Station data from Nordic Arctic Research Program}
\description{Nordic Arctic Research Program (NARP) Long-term variations
  in atmospheric circulation and climate in the Arctic. These data have
  been freely available from a server at met.no, but due to
  re-organisation, they are no longer on the web. Thanks to Ole Einar
  Tveito, I got a copy of the data, which are now inluded in
  clim.pact. These data are accessed theouth \link{getnarp}.}
\usage{data(narp)
       data(narp.meta)
       data(nacd.meta)
       data(nordklim.meta)}
\format{narp is a matrix with 16 columns and narp.meta is a list
  object. The 16th column in narp indicates the type of element:
  code '101', Mean monthly air-temperature; '111', Mean maximum monthly air-temperature;
  '112', Absolute maximum monthly air-temperature; '121', Mean minimum monthly air-temperature;
  '122', Aboslute minimum air-temperature; '401', Mean monthly sea level pressure;
  '601', Mean monthly precipitation sum; '602', Highest monthly 1-day precipitation;
  '701', Mean monthly days with snow cover > 50 percent; '801', Mean monthly cloud cover.}
\examples{
data(narp.meta)
data(nacd.meta)
}
\source{Eirik Forland and Ole Einar Tveito, The Norwegian Meteorological
  Institute; Trausti Jonsson, Icelandic Met Office; Claus Kern-Hansen, Data & Climate
  Division, Danish Meteorological Institute}
\references{
  Ellen Vaarby Laursen (2003) 'DMI monthly climate data, 1873-2002,
  contribution to Nordic Arctic Research Programme (NARP)', DANISH
  METEROLOGICAL INSTITUTE, TECHNICAL REPORT 03-25
  (http://www.dmi.dk/dmi/tr03-25.pdf);
  http://thule.oulu.fi/narp/Projects/a_natural/Forland.htm}
\keyword{datasets}
