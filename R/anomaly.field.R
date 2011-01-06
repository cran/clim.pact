# R.E. Benestad

anomaly.field <- function(x,period=NULL) {
  if ((class(x)[2]!="monthly.field.object") & (class(x)[2]!="field.object") &
      (class(x)[2]!="daily.field.object") & (class(x)[1]!="field")) {
      stop("Need a field.object") }
  nx <- length(x$lon)
  ny <- length(x$lat)
  nt <- length(x$tim)
  i.yy <- is.finite(x$yy)
  if (!is.null(period)) {
    i.yy <- (x$yy >= period[1]) & (x$yy <= period[2])
  }
  dd.rng <- range(x$dd)
  if ( (lower.case(substr(attr(x$tim,"unit"),1,5))=="month") |
       ((dd.rng[2]-dd.rng[1]<4) & (x$mm[2]-x$mm[1]>0)) ) {
    clim <- matrix(x$dat[1,,]*NA,ny,nx)
    for (im in 1:12) {
      it <- mod(1:nt,12)==mod(im,12) & i.yy
      for (j in 1:ny) {
        for (i in 1:nx) {
          clim[j,i] <- mean(x$dat[it,j,i],na.rm=TRUE)
          x$dat[it,j,i] <- x$dat[it,j,i] - clim[j,i]
        }
      }
    }
  } else {
    nt <- sum(i.yy)
    time <- julday(x$mm,x$dd,x$yy) - julday(1,1,1950)
    if (!is.null(x$attributes$daysayear)) daysayear <- x$attributes$daysayear else
                                          daysayear <- 365.25
    x.1<-cos(2*pi*time/daysayear)
    x.2<-sin(2*pi*time/daysayear)
    x.3<-cos(4*pi*time/daysayear)
    x.4<-sin(4*pi*time/daysayear)
    x.5<-cos(6*pi*time/daysayear)
    x.6<-sin(6*pi*time/daysayear)
    dim(x$dat) <- c(nt,ny*nx)
    for (ip in seq(1,ny*nx,by=1)) {
      if (sum(is.finite(x$dat[,ip])) > 0) {
        calibrate <- data.frame(y=x$dat[i.yy,ip],x1=x.1[i.yy],x2=x.2[i.yy],
                                x3=x.3[i.yy],x4=x.4[i.yy],x5=x.5[i.yy],x6=x.6[i.yy])
        ac.fit<-lm(y ~ x1 + x2 + x3 + x4 + x5 + x6, data=calibrate)
        ac <- data.frame(x1=x.1,x2=x.2,x3=x.3,x4=x.4,x5=x.5,x6=x.6)
        x$dat[,ip]<- x$dat[,ip] - predict(ac.fit,newdata=ac)
      } else x$dat[,ip]<- rep(NA,nt)
    }
  dim(x$dat) <- c(nt,ny,nx)
  }

  attr(x$dat,'description') <- 'anomaly'
  invisible(x)  
}


daily2monthly.field <- function(field,min.days.month=20,method="colMeans",na.rm=TRUE) {
  if (lower.case(class(field)[2])!="daily.field.object") stop("daily2monthly.field: needsdaily.field.object ")
  x <- field$dat

  unit <- field$unit
  ele <- field$ele
  ny <- length(field$lat); nx <- length(field$lon); nt <- length(field$tim)
  years <- as.numeric(rownames(table(field$yy))); nyrs <-  length(years)
  months <- 1:12
  val <- rep(NA,nyrs*12*ny*nx)
  dim(val) <- c(nyrs*12,ny*nx)
  dim(x) <- c(nt, ny*nx)
  mm <- matrix(rep(NA,12*nyrs),nyrs,12); yy <- mm
  for (iy in 1:nyrs) {
    for (im in 1:12) {
      this.month.year <-  is.element(field$mm,im) & is.element(field$yy,years[iy])
      if (sum(this.month.year,na.rm=TRUE) > min.days.month) {
        if (!na.rm) val[iy,im] <-
              eval(parse(text=paste(method,"(x[this.month.year,])",sep=""))) else 
            val[(iy-1)*12+im,] <-
              eval(parse(text=paste(method,"(x[this.month.year,],na.rm=TRUE)",sep="")))
            mm[iy,im] <- im;  yy[iy,im] <- years[iy]
      } else {
        print(paste(im,years[iy],"number of days=",sum(this.month.year,na.rm=TRUE)))
        mm[iy,im] <- NA; yy[iy,im] <- NA
      }
    }
  }
  mm <- c(t(mm)); yy <- c(t(yy))
  good <- is.finite(mm) & is.finite(yy)
  val <- val[good,]; mm <- mm[good]; yy <- yy[good]
  plot#(yy); stop()
  nt <- sum(good)
  dim(val) <- c(nt,ny,nx)
  field$attributes$time.unit <- "month"
  res <-  list(dat=val,lon=field$lon,lat=field$lat,tim=1:nt,lev=field$lev,
               v.name=field$v.name,id.x=field$id.x,id.t=rep(field$id.t[1],nt),
               yy=yy,mm=mm,dd=rep(15,nt),n.fld=1,
               id.lon=rep(field$v.name,nx),id.lat=rep(field$v.name,ny),
               attributes=c(field$attributes,'daily2monthly.field'),filename=field$filename)
  class(res) <- c("field","monthly.field.object")
  invisible(res)
}
