plotField <- function(x,lon=NULL,lat=NULL,tim=NULL,mon=NULL,
                       col="black",lty=1,lwd=1,what="ano") {

  if ((class(x)[1]!="field") & (class(x)[1]!="monthly.field.object") &
      (class(x)[1]!="daily.field.object")){
    stop("x must be a 'field' object.")
  }

  if (is.null(lon) & is.null(lat) & is.null(tim)) {
    stop("At least one of lon/lat/tim must be specified for 2D plots.")
  }
  
 cmon<-c('Jan','Feb','Mar','Apr','May','Jun',
         'Jul','Aug','Sep','Oct','Nov','Dec')
  
  if (!is.null(mon)) {
    im <- x$mm== mon
    x$dat <- x$dat[im,,]
    x$yy <- x$yy[im]
    x$mm <- x$mm[im]
    x$dd <- x$dd[im]
    x$id.t <- x$id.t[im]
    date <- cmon[mon]
  }


  # Hovmuller diagrams
  dims <- dim(x$dat)
  nt <- dims[1]
  lon.tim=FALSE
  tim.lat=FALSE
  lon.lat=FALSE
  time.ts=FALSE
  if (!is.null(lon) & is.null(lat) & is.null(tim)) {
    dx <- mean(diff(x$lon),na.rm=TRUE)
    ii <- (x$lon >= min(lon)) & (x$lon < max(lon) + dx)
    if (sum(ii)==1) Z <- x$dat[,,ii] else if (sum(ii)>1) {
      Z <- x$dat[,,1]*0
      for (j in 1:dims[2]) {
        for (i in 1:dims[1]) {
          Z[i,j] <- mean(x$dat[i,j,ii],na.rm=TRUE)
        }
      }
    }
    xlab <- "Time"
    ylab <- "Latitude (deg N)"
    X <- x$yy + (x$mm - 0.5)/12
    Y <- x$lat
    tim.lat=TRUE
    np <- dims[2]
    dim(Z) <- c(dims[1],dims[2])
  }
  if (is.null(lon) & !is.null(lat) & is.null(tim)) {
    dy <- mean(diff(x$lat),na.rm=TRUE)
    ii <- (x$lat >= min(lat)) & (x$lat < max(lat) + dy)
    if (sum(ii)==1) Z <- t(x$dat[,ii,]) else if (sum(ii)>1) {
      Z <- x$dat[,,1]*0
      for (j in 1:dims[1]) {
        for (i in 1:dims[3]) {
          Z[i,j] <- mean(x$dat[j,ii,i],na.rm=TRUE)
        }
      }
    }
    ylab <- "Time"
    xlab <- "Longitude (deg E)"
    Y <- x$yy + (x$mm - 0.5)/12
    X <- x$lon
    lon.tim=TRUE
    np <- dims[3]
    dim(Z) <- c(dims[3],dims[1])
  }
  # Map - call lower level plot function:
  if (is.null(lon) & is.null(lat) & !is.null(tim)) {
    ii <- (x$yy >= min(tim)) & (x$yy < max(tim) + 1)
    if (sum(ii) == 1) {
      l <- seq(1,dims[1],by=1)[ii]
      x$tim <- x$tim[im]
      results <- map.field(x,l=l,what=what,col=col,lty=lty,lwd=lwd)
    } else if (sum(ii) > 1) {
      map <- mean.field(x,tim=x$yy[ii])
      results <- map.map(x,col=col,lty=lty,lwd=lwd)
    }
    lon.lat <- TRUE
  }
  # Time series - call lower level plot function:
  if (!is.null(lon) & !is.null(lat) & is.null(tim)) {
    results <- grd.box.ts(x,lon,lat,what=what,col=col,lty=lty,lwd=lwd)
    time.ts <- TRUE
  }

  if ((tim.lat) | (lon.tim)) {
    if (is.null(mon) & (what=="ano")) {
      if ( (x$mm[2]-x$mm[1]>=1) & (x$dd[2]==x$dd[1]) ) {
      for (im in 1:12) {
        it <- x$mm==im
        clim <- Z*0
        for (ip in 1:np) {
           if (tim.lat) clim[it,ip] <- mean(Z[it,ip],na.rm=TRUE)
           if (lon.tim) clim[ip,it] <- mean(Z[ip,it],na.rm=TRUE)
        }
      }
    } else {
        ac.mod<-matrix(rep(NA,nt*6),nt,6)
        ac.mod[,1]<-cos(2*pi*x$tim/365.25)
        ac.mod[,2]<-sin(2*pi*x$tim/365.25)
        ac.mod[,3]<-cos(4*pi*x$tim/365.25)
        ac.mod[,4]<-sin(4*pi*x$tim/365.25)
        ac.mod[,5]<-cos(6*pi*x$tim/365.25)
        ac.mod[,6]<-sin(6*pi*x$tim/365.25)
        dim(x$dat) <- c(nt,np)
        for (ip in seq(1,np,by=1)) {
          if (tim.lat) vec <- Z[,ip] else vec <- Z[ip,]
          if (sum(is.finite(vec)) > 0) {
            ac.fit<-lm(vec ~ ac.mod)
            if (tim.lat) clim[,ip] <- ac.fit$fit else clim[ip,] <- ac.fit$fit
          } else {
            if (tim.lat) clim[,ip] <- NA else clim[ip,] <- NA
          }
        }
      }

      Z <- switch(lower.case(substr(what,1,3)),
                  "ano"=Z - clim,
                  "cli"=clim,
                  "abs"=Z)
  }

    z.levs <- seq(min(as.vector(Z),na.rm=TRUE),
                  max(as.vector(Z),na.rm=TRUE),length=41)
    my.col <- rgb(c(seq(0,1,length=20),rep(1,21)),
                  c(abs(sin((0:40)*pi/40))),
                  c(c(rep(1,21),seq(1,0,length=20))))
    filled.contour(X,Y,Z,
                   col = my.col,levels=z.levs,
                   main=paste(attributes(x$dat)$"long_name"),
                   sub=date,xlab=xlab,ylab=ylab)

# From filled.contour in base
    mar.orig <- (par.orig <- par(c("mar","las","mfrow")))$mar
    on.exit(par(par.orig))

    w <- (3 + mar.orig[2]) * par('csi') * 2.54
    layout(matrix(c(2, 1), nc=2), widths=c(1, lcm(w)))
    
    par(las = 1)
    mar <- mar.orig
    mar[4] <- 1
    par(mar=mar)
    contour(X,Y,Z,add=TRUE,col=col,lwd=lwd,lty=lty)
    results <- list(Z=Z,x=X,y=Y,xlab=xlab,ylab=ylab,
                    descr=paste(attributes(x$dat)$"long_name"))
  }
  class(results) <- "2D field"
  attr(results,"long_name") <- attr(x$dat,"long_name")
  attr(results,"descr") <- "plotField"
  invisible(results)
}  
