\name{avail.elem}
\alias{avail.elem}
\alias{avail.locs}
\alias{avail.preds}
\alias{avail.eofs}
\alias{avail.ds}
\title{Available elements}
\description{The function \code{avail.preds} searches for available predictors
  and returns a list of file names. \code{avail.elem} returns a list of
  available elements and \code{avail.locs} a list of available locations. These
  functions are support functions for \code{\link{getnacd}} and \code{\link{getnordklim}}.}
\usage{
avail.ds(direc="output")
avail.eofs(direc="data")
avail.preds(direc="data")
avail.elem()
avail.locs(ele) }
\arguments{
    \item{direc}{String containing the data directory.} 
    \item{ele}{Integer code for element in the Nordklim dataset:}
    \tabular{ll}{
      101 \tab mean T(2m) \cr
      111 \tab mean maximum T(2m) \cr
      112 \tab highest maximum T(2m) [Th] \cr
      113 \tab day of Th date Thd \cr
      121 \tab mean minimum T(2m) \cr
      122 \tab lowest minimum T(2m) [Tl] \cr
      123 \tab day of Tl date Tld \cr
      401 \tab mean SLP \cr
      601 \tab monthly accum. precip. \cr
      602 \tab maximum precip. \cr
      701 \tab Number of days with snow cover (> 50\% covered) days dsc \cr
      801 \tab Mean cloud cover \% N  \cr
      911 \tab mean snow depth \cr}
}
\value{
    
    \tabular{ll}{
      avail.preds \tab vector of characters \cr
      avail.locs \tab a list with name, lons, lats,
      country, ident \cr
      avail.elem \tab a list with data.set, ele, name \cr}}
\author{R.E. Benestad}}
\examples{
library(clim.pact)
avail.elem()$name
# [1] "mean T(2m)"                                             
# [2] "mean maximum T(2m)"                                     
# [3] "highest maximum T(2m)"                                  
# [4] "day of Th date Thd"                                     
# [5] "mean minimum T(2m)"                                     
# [6] "lowest minimum T(2m)"                                   
# [7] "day of Tl date Tld"                                     
# [8] "mean SLP"                                               
# [9] "monthly accum. precip."                                 
#[10] "maximum precip."                                        
#[11] "Number of days with snow cover (> 50% covered) days dsc"
#[12] "Mean cloud cover % N"                                   
#[13] "mean snow depth"

# The following assumes that the subdirectory 'data' exists
\dontrun{
avail.locs()$name[avail.locs()$country=="FIN"]
# [1] "HELSINKI"      "TURKU"         "TAMPERE"       "LAPPEENRANTA" 
# [5] "JYVASKYLA"     "KUOPIO"        "KAJAANI"       "OULU"         
# [9] "KUUSAMO"       "SODANKYLA"     "Maarianhamina" "Helsinki"     
#[13] "Turku"         "Huittinen"     "Tampere"       "Hattula"      
#[17] "Heinola"       "Virolahti"     "Lappeenranta"  "Lavia"        
#[21] "Virrat"        "Orivesi"       "Jyvaeskylae"   "Vaasa"        
#[25] "Ylistaro"      "Aehtaeri"      "Kuopio"        "Maaninka"     
#[29] "Joensuu"       "Kestilä"       "Kajaani"       "Oulu"         
#[33] "Yli-Ii"        "Pudasjärvi"    "Kuusamo"       "Sodankylae"

avail.preds()
# [1] "eof.dc.Rdata"                                                 
# [2] "eof.dmc.Rdata"                                                
# [3] "eof.mc2.Rdata"                                                
# [4] "eof.nn.dc.Rdata"                                              
# [5] "eof.nn.dmc.Rdata"
# ...
}
}
\keyword{file}
