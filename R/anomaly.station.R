# R.E. Benestad, met.no, Oslo, Norway 22.05.2002
# rasmus.benestad@met.no
#-------------------------------------------------------------------
# Estimate anomalies

anomaly.station <- function(obs,period=c(1961,1990)) {


cmon<-c("Jan","Feb","Mar","Apr","May","Jun",
        "Jul","Aug","Sep","Oct","Nov","Dec")

if (lower.case(class(obs)[2])=="monthly.station.record") {
  ny <- length(obs$yy)
  value <- t(obs$val)
  if (!is.null(period)) ii <- ((obs$yy>=period[1]) & (obs$yy<=period[2])) else
                        ii <- is.finite(obs$yy)
  for (im in 1:12) {
        value[im,] <- value[im,] - mean(value[im,ii],na.rm=TRUE)
      }
  obs$val <- t(value)
  obs$obs.name  <-  paste(obs$obs.name,"anomaly")
  }  else if (lower.case(class(obs)[2])=="daily.station.record") {

    if (!is.null(attr(obs$tim,"daysayear"))) daysayear <- attr(obs$tim,"daysayear") else
                                             daysayear <- 365.25
    nt <- length(oslo$t2m)
    ac.mod<-matrix(rep(NA,nt*6),nt,6)
    jtime <- 1:nt
    ac.mod[,1]<-cos(2*pi*jtime/daysayear); ac.mod[,2]<-sin(2*pi*jtime/daysayear)
    ac.mod[,3]<-cos(4*pi*jtime/daysayear); ac.mod[,4]<-sin(4*pi*jtime/daysayear)
    ac.mod[,5]<-cos(6*pi*jtime/daysayear); ac.mod[,6]<-sin(6*pi*jtime/daysayear) 
    #print(c(dim(ac.mod),length(obs$t2m)))
    ac.fit <- lm(obs$t2m ~ ac.mod)
    clim <- obs$t2m
    clim[is.finite(obs$t2m)]<-ac.fit$fit
    #print(c(dim(ac.mod),length(obs$t2m), length(clim),NA,sum(!is.finite(obs$t2m))))
    obs$abs.t2m <- obs$t2m
    obs$clim.t2m <- clim
    obs$t2m <- obs$t2m - clim
    ac.fit <- glm(obs$precip ~ ac.mod)
    clim <- obs$precip
    clim[is.finite(obs$precip)]<-ac.fit$fit
    obs$abs.precip <- obs$precip
    obs$clim.precip <- clim
    obs$precip <- obs$precip - clim
  }
  invisible(obs)
}
