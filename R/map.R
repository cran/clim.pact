map <- function(x,y=NULL,col="black",lwd=1,lty=1,sym=TRUE,
                    plot=TRUE,inv.col=FALSE,add=FALSE,las = 1) {

  if (!is.null(y)) map <- x$map - y$map else  map <- x$map

  if (sym) {
    z.levs <- seq(-max(abs(as.vector(map)),na.rm=TRUE),
                   max(abs(as.vector(map)),na.rm=TRUE),length=41)
  } else {
    z.levs <- seq(min(as.vector(map),na.rm=TRUE),
                  max(as.vector(map),na.rm=TRUE),length=41)
  }
  if (plot) {
    my.col <- rgb(c(seq(0,1,length=20),rep(1,21)),
                  c(abs(sin((0:40)*pi/40))),
                  c(c(rep(1,21),seq(1,0,length=20))))
    if (inv.col) my.col <- reverse(my.col)
    if (!add) { filled.contour(x$lon,x$lat,map,
                               col = my.col,levels=z.levs,
                               main=paste(attributes(x$dat)$"long_name"),
                               sub=x$date,xlab="Longitude",ylab="Latitude")
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
    contour(x$lon,x$lat,map,add=TRUE,col=col,lwd=lwd,lty=lty)
    addland()
  }
  results <- list(lon=x$lon,lat=x$lat,map=map,v.name=,x$v.name,
                  tim=NULL,date=NULL,attributes=x$attributes)
  class(results) <- "map"
#  attr(results) <- attr(x)
  attr(results,"descr") <- "Mean values"
  invisible(results) 
}
