# R.E. Benestad

anomaly.field <- function(x,period=NULL) {
  if ((class(x)!="monthly.field.object") & (class(x)!="field.object") &
      (class(x)!="daily.field.object") & (class(x)!="field")) {
      stop("Need a field.object") }
  nx <- length(x$lon)
  ny <- length(x$lat)
  nt <- length(x$tim)
  i.yy <- is.finite(x$yy)
  if (!is.null(period)) {
    i.yy <- (x$yy >= period[1]) & (x$yy <= period[2])
  }
  if ( (x$mm[2]-x$mm[1]>=1) & (x$dd[2]==x$dd[1]) ) {
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
    time <- julian(x$mm,x$dd,x$yy)
    x.1<-cos(2*pi*time/365.25)
    x.2<-sin(2*pi*time/365.25)
    x.3<-cos(4*pi*time/365.25)
    x.4<-sin(4*pi*time/365.25)
    x.5<-cos(6*pi*time/365.25)
    x.6<-sin(6*pi*time/365.25)
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
