# R.E. Benestad

mapField <- function(x,l=NULL,greenwich=TRUE,what="ano",method="nice",
                       col="black",lwd=2,lty=1,add=FALSE,las = 1) {
  if ((class(x)!="monthly.field.object") & (class(x)!="field.object") &
      (class(x)!="daily.field.object") & (class(x)!="field")) {
      stop("Need a field.object") }
  if (is.null(l)) l <- length(x$tim)
  nx <- length(x$lon)
  ny <- length(x$lat)
  nt <- length(x$tim)
  clim <- x$dat[l,,]
  if ( (x$mm[2]-x$mm[1]>=1) & (x$dd[2]==x$dd[1]) ) {
    it <- mod(1:nt,12)==mod(l,12)
    for (j in 1:ny) {
      for (i in 1:nx) {
        clim[j,i] <- mean(x$dat[it,j,i],na.rm=TRUE)
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
      dim(x$dat) <- c(nt,ny*nx)
      dim(clim) <- c(ny*nx)
      for (ip in seq(1,ny*nx,by=1)) {
        if (sum(is.finite(x$dat[,ip])) > 0) {
          ac.fit<-lm(x$dat[,ip] ~ ac.mod)
          clim[ip]<-ac.fit$fit[l]
        } else clim[ip]<- NA
      }
      dim(x$dat) <- c(nt,ny,nx)
      dim(clim) <- c(ny,nx)
    }

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
  z.levs <- seq(-max(abs(as.vector(map)),na.rm=TRUE),
                 max(abs(as.vector(map)),na.rm=TRUE),length=41)
  my.col <- rgb(c(seq(0,1,length=20),rep(1,21)),
                c(abs(sin((0:40)*pi/40))),
                c(c(rep(1,21),seq(1,0,length=20))))
  if ((!add) & (method!="nice")) {
        image(x$lon,x$lat,t(map),levels=z.levs,
        main=paste(attributes(x$dat)$"long_name",descr),
        sub=date,xlab="Longitude",ylab="Latitude")
       } else {
         filled.contour(x$lon,x$lat,t(map),
                        col = my.col,levels=z.levs,
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
  contour(x$lon,x$lat,t(map),add=TRUE,col=col,lwd=lwd,lty=lty)
  addland()
  results <- list(map=t(map),lon=x$lon,lat=x$lat,tim=x$tim,
                  date=date,description=descr,attributes=x$attributes)
  class(results) <- "map"
  attr(results,"long_name") <- attr(x$dat,"long_name")
  attr(results,"descr") <- descr
  invisible(results)

}
