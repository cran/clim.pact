# R.E. Benestad, met.no, Oslo, Norway 22.05.2002
# rasmus.benestad@met.no
#-------------------------------------------------------------------
# Estimate anomalies

anomaly.station <- function(obs,period=c(1961,1990)) {


cmon<-c("Jan","Feb","Mar","Apr","May","Jun",
        "Jul","Aug","Sep","Oct","Nov","Dec")

if (lower.case(class(obs))=="monthly.station.record") {
  ny <- length(obs$yy)
  value <- t(obs$val)
  if (!is.null(period)) ii <- ((obs$yy>=period[1]) & (obs$yy<=period[2])) else
                        ii <- is.finite(obs$yy)
  for (im in 1:12) {
        value[im,] <- value[im,] - mean(value[im,ii],na.rm=TRUE)
      }
  obs$val <- t(value)
  obs$obs.name  <-  paste(obs$obs.name,"anomaly")
  invisible(obs)
  }  
}
