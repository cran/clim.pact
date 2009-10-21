grd.box.ts <- function(x,lon,lat,lev=NULL,what="abs",greenwich=TRUE,mon=NULL,
                       col="grey10",lwd=1,lty=1,pch=".",add=FALSE,
                       filter=NULL,type="l",main=NULL,sub=NULL,xlab=NULL,ylab=NULL,
                       xlim=NULL,ylim=NULL) {

  library(akima)
#  library(ts)
  
  if ((class(x)[1]!="field") & (class(x)[2]!="monthly.field.object") &
      (class(x)[2]!="daily.field.object") ) stop("Need a field.object")

  n.dims <- length(dim(x$dat))
  if ( (n.dims==4) & (is.null(lev)) )  stop("For 4D objects, the level must be given")
  if (greenwich) {
    x$lon[x$lon > 180] <- x$lon[x$lon > 180]-360
    x.srt <- order(x$lon)
    x$lon <- x$lon[x.srt]
    #print(n.dims); print(dim(x$dat)); print(length(x.srt))
    if (n.dims==3) x$dat <- x$dat[,,x.srt] else
    if (n.dims==4) x$dat <- x$dat[,,,x.srt]
  }
  daysayear <- 365.25
  cmon<-c('Jan','Feb','Mar','Apr','May','Jun',
          'Jul','Aug','Sep','Oct','Nov','Dec')
  descr <- "Interpolated value"
  date <- " "
  if (!is.null(mon)) {
    im <- x$mm== mon
    if (n.dims==3) x$dat <- x$dat[im,,] else
    if (n.dims==4) x$dat <- x$dat[im,,,]
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
  n.dims <- length(dim(x$dat))
  if (n.dims==4) {
      if (length(x$lev)>1) {
        dz <- x$lev[2] - x$lev[1]
        z.keep <- (1:length(x$lev))[(x$lev >= lev)][1] 
        x$lev <- x$lev[z.keep]
      } else if (length(x$lev)==1){
        z.keep <- 1
        x$lev <- x$lev[z.keep]
        #print(dim(x$dat))
      }  
  }
  x$lon <- x$lon[x.keep]
  x$lat <- x$lat[y.keep]
  if (n.dims==3) x$dat <- x$dat[,y.keep,x.keep] else
  if (n.dims==4) x$dat <- x$dat[,z.keep,y.keep,x.keep]
  if (sum(!is.finite(x$dat))>0) x$dat[!is.finite(x$dat)] <- 0
  lat.x<-rep(x$lat,length(x$lon))
  lon.x<-sort(rep(x$lon,length(x$lat)))
  nt <- length(x$yy)
  y <- rep(NA,nt)
  for (it in 1:nt) {
    if (n.dims==3) Z.in<-as.matrix(x$dat[it,,]) else
      if (n.dims==4) Z.in<-as.matrix(x$dat[it,z.keep,,])
    Z.out<-interp(lat.x,lon.x,Z.in,lat,lon)
    y[it] <- Z.out$z
  }
#  print("time unit")
  
if (!is.null(attributes(x$tim)$unit)) {
  attr(x$tim,"units") <- attributes(x$tim)$unit
}
  #print(attributes(x$tim)$units)
  #print(attributes(x$tim)$unit)
  #print(summary(y))

  tunit <- attributes(x$tim)$units
  if (!is.null(tunit)) tunit <- lower.case(substr(tunit,1,3)) else {
        tunit <- attributes(x$tim)$units
        if (!is.null(tunit)) tunit <- lower.case(substr(tunit,1,3)) else
                             if (min(diff(x$mm))==1) tunit <- "mon" else
                                                     tunit <- "day"
  }
                            
  if (tunit== "mon") {
    clim <- y
    for (im in 1:12) {
      ii <- mod((1:nt)-1,12)+1 == im
      clim[ii] <- mean(y[ii],na.rm=TRUE)
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

  #print("what?")
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

  #print(summary(ts)); print("plot")
  if (!add) {
  
     plot(x$yy+x$mm/12+x$dd/daysayear,ts,type=type,pch=pch,xlim=xlim,ylim=ylim,
       main=main,sub=sub,xlab=xlab,ylab=ylab,col=col,lwd=lwd,lty=lty)
     points(x$yy+x$mm/12+x$dd/daysayear,ts,pch=pch,col=col)
   } else {
     if (type!='p') lines(x$yy+x$mm/12+x$dd/daysayear,ts,type=type,col=col,lwd=lwd,lty=lty)
     points(x$yy+x$mm/12+x$dd/daysayear,ts,pch=pch,col=col)
   }
  grid()
  #print("plotted")
  
  dd.rng <- range(x$dd)
  if (is.null(attr(x$tim,"units"))) attr(x$tim,"units") <- "unknown"
  if ( (tunit=="mon") |
       ((dd.rng[2]-dd.rng[1]<4) & (x$mm[2]-x$mm[1]>0)) ) {
#    print("Monthly")
    results <- station.obj(ts,yy=x$yy,mm=x$mm,obs.name=x$v.name,
                           unit=x$attributes$unit,ele=NA,
                           station=NA,lat=round(lat,4),lon=round(lon,4),alt=NA,
                           location="interpolated",wmo.no=NA,
                           start=min(x$yy),yy0=attr(x$tim,"time_origin"),country=NA,
                           ref="grd.box.ts.R (clim.pact)")
  } else {
    attr(x$tim,"daysayear") <- daysayear
    results <- station.obj.dm(t2m=ts,precip=rep(NA,length(ts)),
                              dd=x$dd,mm=x$mm,yy=x$yy,
                              obs.name=x$v.name,unit=x$attributes$unit,ele=NA,
                              station=NA,lat=round(lat,4),lon=round(lon,4),alt=NA,
                              location="interpolated",wmo.no=NA,
                              start=min(x$yy),yy0=attr(x$tim,"time_origin"),country=NA,
                              ref="grd.box.ts.R (clim.pact)")
  }

#  print("exit grd.box.ts()")
  invisible(results)
}


