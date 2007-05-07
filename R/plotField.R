plotField <- function(x,lon=NULL,lat=NULL,tim=NULL,mon=NULL,val.rng=NULL,
                      col="black",col.coast="grey",lty=1,lwd=1,what="ano",
                      type="l",pch=26,my.col=NULL,add=FALSE,
                      main=NULL,sub=NULL,xlab=NULL,ylab=NULL,
                      xlim=NULL,ylim=NULL) {

  if ((class(x)[1]!="field") & (class(x)[1]!="monthly.field.object") &
      (class(x)[1]!="daily.field.object")){
    stop("x must be a 'field' object.")
  }

  if (is.null(lon) & is.null(lat) & is.null(tim)) {
    stop("At least one of lon/lat/tim must be specified for 1D or 2D plots.")
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
  ind <- 1:length(x$tim)

  # Hovmuller diagrams
  dims <- dim(x$dat)
  nt <- dims[1]
  lon.tim=FALSE
  tim.lat=FALSE
  lon.lat=FALSE
  time.ts=FALSE

  if (!is.null(lon) & is.null(lat) & is.null(tim)) {
    #print("Time-lat")
    dx <- mean(diff(x$lon),na.rm=TRUE)
    ii <- (x$lon >= min(lon)) & (x$lon < max(lon) + dx)
    if (sum(ii)==1) Z <- x$dat[,,ii] else if (sum(ii)>1) {
      print("(sum(ii)>1)")
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
  } else {
    # print("did nothing!")
    
  }

  if (is.null(lon) & !is.null(lat) & is.null(tim)) {
    #print("Time-lon")
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
  if ( (is.null(lon) & is.null(lat) & !is.null(tim)) |
       (is.null(lon) & length(lat)==2 & !is.null(tim)) |
       (length(lon)==2 & length(lat)==2 & !is.null(tim)) |
       (length(lon)==2 & is.null(lat) & !is.null(tim)) ){
    print("Map")
    ii <- (ind >= min(tim)) & (ind < max(tim) + 1)
    if ( (length(lon)==2) & is.null(xlim) ) xlim <- lon
    if ( (length(lat)==2) & is.null(ylim) ) ylim <- lat
    if (sum(ii) == 0) {
      ii <- rep(FALSE,length(x$tim))
      ii[tim] <- TRUE
    }
    if (sum(ii) == 1) {
      l <- seq(1,dims[1],by=1)[ii]
#      print(paste("l=",l,"date=",x$yy[ii],x$mm[ii],x$dd[ii]))
      results <- mapField(x,l=l,what=what,col=col,col.coast=col.coast,lty=lty,
                          lwd=lwd,val.rng=val.rng,xlim=xlim,ylim=ylim)
    } else if (sum(ii) > 1) {
      map <- meanField(x,xlim=lon,ylim=lat,t.rng=range(x$yy[ii]))
      results <- mapField(x,col=col,col.coast=col.coast,lty=lty,lwd=lwd,val.rng=val.rng)
    }       
    lon.lat <- TRUE
    invisible(results)
    return()
  } 
  # Time series - call lower level plot function:
  if (!is.null(lon) & !is.null(lat)) {
    #print("plotField: Time-series")
    results <- grd.box.ts(x,lon,lat,what=what,col=col,
                          lty=lty,lwd=lwd,pch=pch,type=type,add=add,
                          main=main,sub=sub,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
    Z <- results$t2m
    time.ts <- TRUE
#    print("Now, try to exit this...")
  } 

  if ((tim.lat) | (lon.tim)) {
    print("Hovmuller diagrams")
    clim <- Z*0
    if (is.null(mon) & (what=="ano")) {
      if ( (x$mm[2]-x$mm[1]>=1) & (x$dd[2]==x$dd[1]) ) {
      for (im in 1:12) {
        it <- x$mm==im
           for (ip in 1:np) {
           if (tim.lat) clim[it,ip] <- mean(Z[it,ip],na.rm=TRUE)
           if (lon.tim) clim[ip,it] <- mean(Z[ip,it],na.rm=TRUE)
        }
      }
    } else if (lon.lat) {
      #print("Longitude-latitude map")
        ac.mod<-matrix(rep(NA,nt*6),nt,6)
        ac.mod[,1]<-cos(2*pi*x$tim/365.25); ac.mod[,2]<-sin(2*pi*x$tim/365.25)
        ac.mod[,3]<-cos(4*pi*x$tim/365.25); ac.mod[,4]<-sin(4*pi*x$tim/365.25)
        ac.mod[,5]<-cos(6*pi*x$tim/365.25); ac.mod[,6]<-sin(6*pi*x$tim/365.25)
        dim(x$dat) <- c(nt,np)
        ac.fit<-lm(x$dat ~ ac.mod)
        clim <- ac.fit$fit

# slow ...
#        for (ip in seq(1,np,by=1)) {
#          if (tim.lat) vec <- Z[,ip] else vec <- Z[ip,]
#          if (sum(is.finite(vec)) > 0) {
#            ac.fit<-lm(vec ~ ac.mod)
#            if (tim.lat) clim[,ip] <- ac.fit$fit else clim[ip,] <- ac.fit$fit
#          } else {
#            if (tim.lat) clim[,ip] <- NA else clim[ip,] <- NA
#          }
#        }
#--------------------------------------------------------
      }
    }
  }

  if (!time.ts) {
    #print(paste("2D-plots",what))
    Z <- switch(lower.case(substr(what,1,3)),
                              "ano"=Z - clim,
                              "cli"=clim,
                              "abs"=Z)
 
   if (is.null(main)) main <-  paste(attributes(x$dat)$"long_name")
   if (is.null(sub)) sub <- date
    
    if (is.null(val.rng)) {
      z.levs <- seq(min(abs(as.vector(Z)),na.rm=TRUE),
                    max(abs(as.vector(Z)),na.rm=TRUE),length=21)
    } else z.levs <- seq(val.rng[1],val.rng[2],length=21)
    
    if (is.null(my.col)) my.col <- rgb(c(seq(0,1,length=10),rep(1,11)),
                                       c(abs(sin((0:20)*pi/20))),
                                       c(c(rep(1,11),seq(1,0,length=10))))
    if (is.null(xlim)) xlim <- range(X,na.rm=TRUE)
    if (is.null(ylim)) ylim <- range(Y,na.rm=TRUE)
    filled.contour(X,Y,Z,xlim=xlim,ylim=ylim,
                   col = my.col,levels=z.levs,
                   main=main,sub=sub,xlab=xlab,ylab=ylab)

# From filled.contour in base
    mar.orig <- (par.orig <- par(c("mar","las","mfrow")))$mar
    on.exit(par(par.orig))

    w <- (3 + mar.orig[2]) * par('csi') * 2.54
    layout(matrix(c(2, 1), nc=2), widths=c(1, lcm(w)))
    
    par(las = 1)
    mar <- mar.orig
    mar[4] <- 1
    par(mar=mar)
#    contour(X,Y,Z,add=TRUE,col=col,lwd=lwd,lty=lty,levels=z.levs)
    results <- list(Z=Z,x=X,y=Y,xlab=xlab,ylab=ylab,
                    descr=paste(attributes(x$dat)$"long_name"))

#  print("Finished")
  class(results) <- "2D field"
  attr(results,"long_name") <- attr(x$dat,"long_name")
  attr(results,"descr") <- "plotField"
  }
invisible(results)
}

