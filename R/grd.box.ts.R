grd.box.ts <- function(x,lon,lat,what="abs",greenwich=TRUE,mon=NULL,
                       col="grey10",lwd=1,lty=1,pch=26,add=FALSE,
                       filter=NULL,type="s") {

  library(akima)
  library(ts)
  
  if ((class(x)[1]!="field") & (class(x)!="monthly.field.object") &
      (class(x)!="daily.field.object") ) stop("Need a field.object")
  
  if (greenwich) {
    x$lon[x$lon > 180] <- x$lon[x$lon > 180]-360
    x.srt <- order(x$lon)
    x$lon <- x$lon[x.srt]
    x$dat <- x$dat[,,x.srt]
  }

  cmon<-c('Jan','Feb','Mar','Apr','May','Jun',
          'Jul','Aug','Sep','Oct','Nov','Dec')
  descr <- "Interpolated value"
  date <- " "
  if (!is.null(mon)) {
    im <- x$mm== mon
    x$dat <- x$dat[im,,]
    x$yy <- x$yy[im]
    x$mm <- x$mm[im]
    x$dd <- x$dd[im]
    x$id.t <- x$id.t[im]
    date <- cmon[mon]
  }
  dx <- x$lon[2] - x$lon[1]
  dy <- x$lat[2] - x$lat[1]
  x.keep <- (x$lon - 3*dx <= lon) & (x$lon + 3*dx >= lon)
  y.keep <- (x$lat - 3*dy <= lat) & (x$lat + 3*dx >= lat)
  x$lon <- x$lon[x.keep]
  x$lat <- x$lat[y.keep]
  x$dat <- x$dat[,y.keep,x.keep]
  lat.x<-rep(x$lat,length(x$lon))
  lon.x<-sort(rep(x$lon,length(x$lat)))
  nt <- length(x$yy)
  y <- rep(NA,nt)
  for (it in 1:nt) {
    Z.in<-as.matrix(x$dat[it,,])
    Z.out<-interp(lat.x,lon.x,Z.in,lat,lon)
    y[it] <- Z.out$z
  }

#  print("time unit")
  
if (!is.null(attributes(x$tim)$unit)) {
  attr(x$tim,"units") <- attributes(x$tim)$unit
}
#  print(attributes(x$tim)$units)
#  print(attributes(x$tim)$unit)

  if (lower.case(substr(attributes(x$tim)$units,1,5))== "month") {
    clim <- y
    for (im in 1:12) {
      ii <- mod((1:nt)-1,12)+1 == im
      clim[ii] <- mean(y[ii],na.rm=T)
    }
  } else {
    ac.mod<-matrix(rep(NA,nt*6),nt,6)
    if (substr(lower.case(attributes(x$tim)$units),1,3)=="day") jtime <- x$tim
    if (substr(lower.case(attributes(x$tim)$units),1,4)=="hour")  jtime <- x$tim/24
    ac.mod[,1]<-cos(2*pi*jtime/365.25); ac.mod[,2]<-sin(2*pi*jtime/365.25)
    ac.mod[,3]<-cos(4*pi*jtime/365.25); ac.mod[,4]<-sin(4*pi*jtime/365.25)
    ac.mod[,5]<-cos(6*pi*jtime/365.25); ac.mod[,6]<-sin(6*pi*jtime/365.25)
    ac.fit<-lm(y ~ ac.mod); clim <- ac.fit$fit
  } 

#  print("what?")
  ts <- switch(lower.case(substr(what,1,3)),
                "ano"=y - clim,
                "cli"=clim,
                "abs"=y)
  descr <- switch(lower.case(substr(what,1,3)),
                "ano"="anomaly",
                "cli"="climatological",
                "abs"="absolute value")

  if (!is.null(filter)) ts <- filter(ts,filter)
  if (!add) {

     plot(x$yy+(x$mm-0.5)/12,ts,type=type,pch=pch,
       main=x$v.name,
          sub=paste("Interpolated at ",lon,"E, ",lat,"N ",date,sep=""),
       xlab="Time",ylab=attributes(x$dat)$unit,col=col,lwd=lwd,lty=lty)

     points(x$yy+(x$mm-0.5)/12,ts,pch=pch,col=col)
   } else {
     if (type!='p') lines(x$yy+(x$mm-0.5)/12,ts,type=type,col=col,lwd=lwd,lty=lty)
     points(x$yy+(x$mm-0.5)/12,ts,pch=pch,col=col)
   }
  grid()
#  print("plotted")
  
  dd.rng <- range(x$dd)
  if (is.null(attr(x$tim,"units"))) attr(x$tim,"units") <- "unknown"
  if ( (lower.case(substr(attr(x$tim,"units"),1,5))=="month") |
       ((dd.rng[2]-dd.rng[1]<4) & (x$mm[2]-x$mm[1]>0)) ) {
#    print("Monthly")
    results <- station.obj(ts,yy=x$yy,obs.name=x$v.name,unit=attr(x$dat,"unit"),
                           ele=NA,mm=x$mm,
                           station=NA,lat=round(lat,4),lon==round(lon,4),alt=NA,
                           location="interpolated",wmo.no=NA,
                           start=min(x$yy),yy0=attr(x$tim,"time_origin"),country=NA,
                           ref="grd.box.ts.R (clim.pact)")
  } else {
    results <- station.obj.dm(t2m=ts,precip=rep(NA,length(ts)),
                              x$dd,x$mm,x$yy,
                              obs.name=x$v.name,unit=attr(x$dat,"unit"),ele=NA,
                              station=NA,lat=round(lat,4),lon=round(lon,4),alt=NA,
                              location="interpolated",wmo.no=NA,
                              start=min(x$yy),yy0=attr(x$tim,"time_origin"),country=NA,
                              ref="grd.box.ts.R (clim.pact)")
  }

#  print("exit grd.box.ts()")
  invisible(results)
}


