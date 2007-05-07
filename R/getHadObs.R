getHadObs <- function(what="CET") {
  url <- switch(upper.case(what),
                "CET"="http://hadobs.metoffice.com/hadcet/cetml1659on.dat",
                "EWP"="http://hadobs.metoffice.com/hadukp/data/monthly/HadEWP_monthly_qc.txt")
  cnames <- c("yy","JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC","ANNUAL")
  skip=switch(upper.case(what),"CET"=7,"EWP"=4)
  obs <- as.list(read.table(file=url,skip=7,col.names=cnames,header=FALSE,as.is=TRUE))
  obs$val <- cbind(obs$JAN,obs$FEB,obs$MAR,obs$APR,obs$MAY,obs$JUN,
                   obs$JUL,obs$AUG,obs$SEP,obs$OCT,obs$NOV,obs$DEC)
  obs$val[obs$val< -99] <- NA
  obs$location <- what
  obs$obs.name <- switch(upper.case(what),"CET"="MONTHLY MEAN CENTRAL ENGLAND TEMPERATURE",
                     "EWP"="Monthly England & Wales precipitation")
  obs$reference <- switch(upper.case(what),
  "CET"="MANLEY (Q.J.R.METEOROL.SOC., 1974), PARKER ET AL. (INT.J.CLIM., 1992), PARKER AND HORTON (INT.J.CLIM., 2005)",
  "EWP"="Wigley & Jones (J.Climatol.,1987), Gregory et al. (Int.J.Clim.,1991),Jones & Conway (Int.J.Climatol.,1997), Alexander & Jones (ASL,2001)")
  obs$ele <- switch(upper.case(what),"CET"=101,"EWP"=601)
  obs$lon <- switch(upper.case(what),"CET"=-0.7333,"EWP"=)
  obs$lat <- switch(upper.case(what),"CET"=51.77,"EWP"=)
  obs$alt <- switch(upper.case(what),"CET"=60,"EWP"=NA)
  obs$unit <- switch(upper.case(what),"CET"="DEGREES C","EWP"="mm total")
  obs$found <- TRUE
  obs$wmo.no <- NA
  obs$start <- min(obs$yy)
  obs$yy0 <- min(obs$yy)
  obs$country <- "U.K."
  class(obs) <- c("station","monthly.station.record")
  invisible(obs)
}
