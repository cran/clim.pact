grd.box.ts <- function(x,lon,lat,what="abs",greenwich=TRUE,mon=NULL,
                       col="grey10",lwd=1,lty=1,pch=26,add=FALSE,
                       filter=NULL,type="l",main=NULL,sub=NULL,xlab=NULL,ylab=NULL,
                       xlim=NULL,ylim=NULL) {

  library(akima)
#  library(ts)
  
  if ((class(x)[1]!="field") & (class(x)[2]!="monthly.field.object") &
      (class(x)[2]!="daily.field.object") ) stop("Need a field.object")
  
  if (greenwich) {
    x$lon[x$lon > 180] <- x$lon[x$lon > 180]-360
    x.srt <- order(x$lon)
    x$lon <- x$lon[x.srt]
    x$dat <- x$dat[,,x.srt]
  }
  daysayear <- 365.25
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
  if (sum(!is.finite(x$dat))>0) x$dat[!is.finite(x$dat)] <- 0
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
  #print(attributes(x$tim)$units)
  #print(attributes(x$tim)$unit)

  tunit <- attributes(x$tim)$units
  if (!is.null(tunit)) tunit <- lower.case(substr(tunit,1,3)) else
                       tunit <- "mon"
                       
  if (tunit== "mon") {
    clim <- y
    for (im in 1:12) {
      ii <- mod((1:nt)-1,12)+1 == im
      clim[ii] <- mean(y[ii],na.rm=T)
    }
  } else {
    ac.mod<-matrix(rep(NA,nt*6),nt,6)
    if (tunit=="day") jtime <- x$tim
    if (tunit=="hou")  jtime <- x$tim/24
    if (!is.null(x$attributes$daysayear)) daysayear <- x$attributes$daysayear else
                                          daysayear <- 365.25
    ac.mod[,1]<-cos(2*pi*jtime/daysayear); ac.mod[,2]<-sin(2*pi*jtime/daysayear)
    ac.mod[,3]<-cos(4*pi*jtime/daysayear); ac.mod[,4]<-sin(4*pi*jtime/daysayear)
    ac.mod[,5]<-cos(6*pi*jtime/daysayear); ac.mod[,6]<-sin(6*pi*jtime/daysayear)
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
  if (is.null(main)) main <-  x$v.name
  if (is.null(sub)) sub <- paste("Interpolated at ",lon,"E, ",lat,"N ",date,sep="")
  if (is.null(xlab)) xlab <- "Time"
  if (is.null(ylab)) ylab <- attributes(x$dat)$unit

  if (!add) {

     plot(x$yy+x$mm/12+x$dd/daysayear,ts,type=type,pch=pch,xlim=xlim,ylim=ylim,
       main=main,sub=sub,xlab=xlab,ylab=ylab,col=col,lwd=lwd,lty=lty)
     points(x$yy+x$mm/12+x$dd/daysayear,ts,pch=pch,col=col)
   } else {
     if (type!='p') lines(x$yy+x$mm/12+x$dd/daysayear,ts,type=type,col=col,lwd=lwd,lty=lty)
     points(x$yy+x$mm/12+x$dd/daysayear,ts,pch=pch,col=col)
   }
  grid()
#  print("plotted")
  
  dd.rng <- range(x$dd)
  if (is.null(attr(x$tim,"units"))) attr(x$tim,"units") <- "unknown"
  if ( (tunit=="mon") |
       ((dd.rng[2]-dd.rng[1]<4) & (x$mm[2]-x$mm[1]>0)) ) {
#    print("Monthly")
    results <- station.obj(ts,yy=x$yy,obs.name=x$v.name,unit=attr(x$dat,"unit"),
                           ele=NA,mm=x$mm,
                           station=NA,lat=round(lat,4),lon==round(lon,4),alt=NA,
                           location="interpolated",wmo.no=NA,
                           start=min(x$yy),yy0=attr(x$tim,"time_origin"),country=NA,
                           ref="grd.box.ts.R (clim.pact)")
  } else {
    attr(x$tim,"daysayear") <- daysayear
    results <- station.obj.dm(t2m=ts,precip=rep(NA,length(ts)),
                              dd=x$dd,mm=x$mm,yy=x$yy,
                              obs.name=x$v.name,unit=attr(x$dat,"unit"),ele=NA,
                              station=NA,lat=round(lat,4),lon=round(lon,4),alt=NA,
                              location="interpolated",wmo.no=NA,
                              start=min(x$yy),yy0=attr(x$tim,"time_origin"),country=NA,
                              ref="grd.box.ts.R (clim.pact)")
  }

#  print("exit grd.box.ts()")
  invisible(results)
}


