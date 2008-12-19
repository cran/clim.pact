# R.E. Benestad

mapField <- function(x,l=NULL,greenwich=TRUE,plot=TRUE,
                     what="ano",method="nice",val.rng=NULL,
                     col="black",col.coast="grey",lwd=2,lty=1,
                     add=FALSE,las = 1,levels=NULL,xlim=NULL,ylim=NULL,newFig=TRUE) {
  if ((class(x)[2]!="monthly.field.object") & (class(x)[2]!="field.object") &
      (class(x)[2]!="daily.field.object") & (class(x)[1]!="field")) {
      stop("Need a field.object") }
  if (is.null(l)) l <- length(x$tim)
  print("mapField here 1")
  nx <- length(x$lon); ny <- length(x$lat); nt <- length(x$tim)
  if (is.character(l)) {
    ldate <- datestr2num(l)
    datematch <- is.element(x$yy*10000+x$mm*100+x$dd,ldate[1]*10000+ldate[2]*100+ldate[3])
    l <- (1:length(x$tim))[datematch]
    print(c(ldate,NA,sum(datematch),range(x$yy),NA,l,NA,dim(x$dat)))
    #print(rbind(x$yy[30:40],x$mm[30:40],x$dd[30:40]))
  }
  clim <- x$dat[l,,]
  dd.rng <- range(x$dd,na.rm=TRUE)
        if (is.null(attr(x$tim,"units"))) attr(x$tim,"units") <- x$attributes$time.unit
  if ( (lower.case(substr(attr(x$tim,"units"),1,5))=="month") |
       ((dd.rng[2]-dd.rng[1]<4) & (x$mm[2]-x$mm[1]>0)) ) {
    it <- mod(1:nt,12)==mod(l,12)
    dim(clim) <- c(ny*nx)
    print(c(nt,ny,nx,NA,dim(x$dat)))
    dim(x$dat) <- c(nt,ny*nx)
    clim <- colMeans(x$dat)
    dim(clim) <- c(ny,nx)
    dim(x$dat) <- c(nt,ny,nx)
# slow...
#    for (j in 1:ny) {
#      for (i in 1:nx) {
#        clim[j,i] <- mean(x$dat[it,j,i],na.rm=TRUE)
#      }
#    }
  } else {
      if (!is.null(x$attributes$daysayear)) daysayear <- x$attributes$daysayear else
                                          daysayear <- 365.25

      ac.mod<-matrix(rep(NA,nt*6),nt,6)
      if (substr(lower.case(attributes(x$tim)$units),1,3)=="day") jtime <- x$tim
      if (substr(lower.case(attributes(x$tim)$units),1,4)=="hour")  jtime <- x$tim/24
      ac.mod[,1]<-cos(2*pi*jtime/daysayear); ac.mod[,2]<-sin(2*pi*jtime/daysayear)
      ac.mod[,3]<-cos(4*pi*jtime/daysayear); ac.mod[,4]<-sin(4*pi*jtime/daysayear)
      ac.mod[,5]<-cos(6*pi*jtime/daysayear); ac.mod[,6]<-sin(6*pi*jtime/daysayear)               
      dim(x$dat) <- c(nt,ny*nx)
      dim(clim) <- c(ny*nx)
      ac.fit <- lm(x$dat ~ ac.mod)
      clim<-ac.fit$fit[l,]
# slow
#      for (ip in seq(1,ny*nx,by=1)) {
#        if (sum(is.finite(x$dat[,ip])) > 0) {
#          ac.fit<-lm(x$dat[,ip] ~ ac.mod)
#          clim[ip]<-ac.fit$fit[l]
#        } else clim[ip]<- NA
#      }
      dim(x$dat) <- c(nt,ny,nx)
      dim(clim) <- c(ny,nx)
    }

  print("mapField here 2")
  cmon<-c('Jan','Feb','Mar','Apr','May','Jun',
          'Jul','Aug','Sep','Oct','Nov','Dec')
  
  if ( (x$mm[2]-x$mm[1]>=1) & (x$dd[2]==x$dd[1]) ) {
    date <- switch(lower.case(substr(what,1,3)),
                   "ano"=paste(cmon[x$mm[l]],x$yy[l]),
                   "cli"=cmon[x$mm[l]],
                   "abs"=paste(cmon[x$mm[l]],x$yy[l]))
  } else {
    date <- switch(lower.case(substr(what,1,3)),
                   "ano"=paste(x$dd[l],cmon[x$mm[l]],x$yy[l]),
                   "cli"=paste(x$dd[l],cmon[x$mm[l]]),
                   "abs"=paste(x$dd[l],cmon[x$mm[l]],x$yy[l]))
  }
  if (greenwich) {
    x$lon[x$lon > 180] <- x$lon[x$lon > 180]-360
    x.srt <- order(x$lon)
    x$lon <- x$lon[x.srt]
    clim <- clim[,x.srt]
    x$dat <- x$dat[,,x.srt]
  }
  anom <- x$dat[l,,]-clim
  map <- switch(lower.case(substr(what,1,3)),
                "ano"=anom,
                "cli"=clim,
                "abs"=x$dat[l,,])
  descr <- switch(lower.case(substr(what,1,3)),
                "ano"="anomaly",
                "cli"="climatological",
                "abs"="absolute value")

  if (plot) {

    if (is.null(val.rng)) {
      print("Mapfield: set range")
      nn <- floor(-max(abs(as.vector(map[is.finite(map)]))))
      xx <- ceiling(max(abs(as.vector(map[is.finite(map)]))))
      nl <- xx-nn
     while (nl > 20) {
        nl <- nl/10
      }
      while (nl < 5) {
        nl <- nl*2
      }
      scl <- 10^floor(log(max(abs(as.vector(map[is.finite(map)]))))/log(10))

      
      print(paste("Scaling is",scl," and range of values is",min(as.vector(map[is.finite(map)])),
                  "-",max(as.vector(map[is.finite(map)])), "[PS. won't plot magnitudes less than 0.001]"))
      if (is.null(levels)) z.levs <- round(seq(nn,xx,length=nl)/scl,2)*scl else {
                           z.levs <- levels; nl <- length(z.levs) }
    print(z.levs)
      my.col <- rgb(c(seq(0,1,length=floor(nl/2)),rep(1,ceiling(nl/2))),
                    c(abs(sin((0:(nl-1))*pi/(nl-1)))),
                    c(c(rep(1,ceiling(nl/2)),seq(1,0,length=floor(nl/2)))))
    print(nl)
    } else {
      z.levs <- seq(val.rng[1],val.rng[2],length=41)
      my.col <- rgb(c(seq(0,1,length=20),rep(1,21)),
                    c(abs(sin((0:40)*pi/40))),
                    c(c(rep(1,21),seq(1,0,length=20))))
    }
    if ((!add) & (method!="nice")) {
        if (newFig) newFig()
        image(x$lon,x$lat,t(round(map,3)),levels=seq(nn,xx,length=101),
        main=paste(attributes(x$dat)$"long_name",descr),
        sub=date,xlab="Longitude",ylab="Latitude")
       } else if (!add) {
         if (newFig) newFig()
         if (is.null(xlim)) xlim <- range(x$lon,na.rm=TRUE)
         if (is.null(ylim)) ylim <- range(x$lat,na.rm=TRUE)

         filled.contour(x$lon,x$lat,t(round(map,3)),
                        col = my.col,levels=z.levs,xlim=xlim,ylim=ylim,
                        main=paste(attributes(x$dat)$"long_name",descr),
                        sub=date,xlab="Longitude",ylab="Latitude")
       }

# From filled.contour in base
    mar.orig <- (par.orig <- par(c("mar","las","mfrow")))$mar
    on.exit(par(par.orig))

    w <- (3 + mar.orig[2]) * par('csi') * 2.54
    layout(matrix(c(2, 1), nc=2), widths=c(1, lcm(w)))
    
    par(las = las)
    mar <- mar.orig
    mar[4] <- 1
    par(mar=mar)
    contour(x$lon,x$lat,t(round(map,3)),add=TRUE,col=col,lwd=lwd,lty=lty,levels=z.levs)
    addland(col=col.coast)
  }   # plot
  
  results <- list(map=t(round(map,3)),lon=x$lon,lat=x$lat,tim=x$tim[l],
                  date=date,description=descr,attributes=x$attributes)
  class(results) <- "map"
  attr(results,"long_name") <- attr(x$dat,"long_name")
  attr(results,"descr") <- descr
  invisible(results)
}
