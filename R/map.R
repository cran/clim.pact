map <- function(x,y=NULL,col="black",lwd=1,lty=1,sym=TRUE,
                    plot=TRUE,inv.col=FALSE,add=FALSE,las = 1,
                    levels=NULL,main=NULL,sub=NULL,xlim=NULL,ylim=NULL) {
  library(akima)

  if (!is.null(y)) {
    nx.1 <- dim(x$map)[1]; ny.1 <- dim(x$map)[2]
    nx.2 <- dim(y$map)[1]; ny.2 <- dim(y$map)[2]
    if ((nx.1==nx.2) & (ny.1==ny.2) & (x$lon[1]=y$lon[1]) & (x$lat[1]=y$lat[1])) {
      map <- x$map - y$map
    } else {
# interpolate:
      lat.y<-sort(rep(y$lat,length(y$lon)))
      lon.y<-rep(y$lon,length(y$lat))
#      image(y$lon,y$lat,y$map,levels=seq(-30,30,by=2),main="map: test"); addland()
      Z.out<-interp(lon.y,lat.y,y$map,x$lon,x$lat)$z
#      contour(x$lon,x$lat,Z.out,levels=seq(-30,30,by=2),add=TRUE,lwd=2)
      map <- x$map - Z.out
#      contour(x$lon,x$lat,map,col="grey",lty=2,add=TRUE); x11()
    }
  }  else  map <- x$map

  if (sym) {
    z.levs <- seq(-max(abs(as.vector(map)),na.rm=TRUE),
                   max(abs(as.vector(map)),na.rm=TRUE),length=41)
  } else {
    z.levs <- seq(min(as.vector(map),na.rm=TRUE),
                  max(as.vector(map),na.rm=TRUE),length=41)
  }
  nc1 <- 20; nc2 <- 21
  if (!is.null(levels)) {
    z.levs <- levels
    nc1 <- floor(length(levels)/2)
    nc2 <- ceiling(length(levels)/2)
  }
  if (is.null(main)) main <- paste(attributes(x$dat)$"long_name")

  if (is.null(x$tim)) if (length(x$tim) > 1) {
    date <- x$tim[1]
    for (i in 2:length(x$tim)) date <- paste(date,"-",x$tim[i],sep="")
    x$tim <- date
  }

  if (is.null(sub)) {
  if (length(x$date)==1) {
    if (is.null(x$tim)) sub <- x$date else
                        sub <- paste(x$date,": ",x$tim,sep="")
  } else {
    if (is.null(x$tim)) sub <- paste(x$date[1],"-", x$date[length(x$date)]) else
                        sub <- paste(x$date[1]," - ", x$date[length(x$date)],": ",x$tim,sep="")
  }
  }

  if (is.null(xlim)) xlim <- range(x$lon[is.finite(x$lon)]) 
  if (is.null(ylim)) ylim <- range(x$lat[is.finite(x$lat)]) 
  if (plot) {
    my.col <- rgb(c(seq(0,1,length=nc1),rep(1,nc2)),
                  c(abs(sin((0:(length(z.levs)-1))*pi/(length(z.levs)-1)))),
                  c(c(rep(1,nc2),seq(1,0,length=nc1))))
    if (inv.col) my.col <- reverse(my.col)
    if (!add) { filled.contour(x$lon,x$lat,map,xlim=xlim,ylim=ylim,
                               col = my.col,levels=z.levs,
                               main=main,sub=sub,xlab="Longitude",ylab="Latitude")
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
    ci <- 10^ceiling(log(max(abs(z.levs)))/log(10))
    cis <- seq(floor(min(z.levs)),ceiling(max(z.levs)),by=0.05*ci)
    contour(x$lon,x$lat,map,add=TRUE,col=col,lwd=lwd,lty=lty,levels=cis)
    addland()
  }
  results <- list(lon=x$lon,lat=x$lat,map=map,v.name=,x$v.name,
                  tim=x$tim,date=x$date,attributes=x$attributes)
  class(results) <- "map"
#  attr(results) <- attr(x)
  attr(results,"descr") <- "Mean values"
  invisible(results) 
}
