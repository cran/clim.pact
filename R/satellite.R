satellite <- function(map.obj,col="black",lwd=2,lty=1,add=FALSE,
                      las = 1,lon.0=0,lat.0=90,
                      ni=100,nj=100, n.nearest=4,max.dist=3) {
  if (class(map.obj)!="map")  stop("Need a map object (map.field)")

  if (is.null(lon.0)) lon.0 <- mean(map.obj$lon)
  if (is.null(lat.0)) lat.0 <- mean(map.obj$lat)
  nx <- length(map.obj$lon)
  ny <- length(map.obj$lat)
  np <- nx*ny
  lats <- rep(map.obj$lat,nx)
  lons <- sort(rep(map.obj$lon,ny))
  map <- t(map.obj$map)
  dim(map) <- c(np,1)
  keep <- (lats-lat.0 >= -90) & (lats-lat.0 <= 90) &
          (lons-lon.0 >= -90) & (lons-lon.0 <= 90)
  map <- map[keep]
  lons <- lons[keep]
  lats <- lats[keep]
  np <- length(map)
  x <- rep(NA,np)
  y <- rep(NA,np)
  r <- sin( pi*(lats-lat.0)/180 )
  x <- r*sin(-pi*(lons-lon.0)/180)
  y <- r*cos(-pi*(lons-lon.0)/180)
#  x.grd <- seq(min(x),max(x),length=ni)
#  y.grd <- seq(min(y),max(y),length=nj)
  x.grd <- seq(-1,1,length=ni)
  y.grd <- seq(-1,1,length=nj)
  dx <- sqrt(mean(diff(x.grd))^2 + mean(diff(y.grd))^2)
  map.xy <- matrix(rep(NA,ni*nj),ni,nj)
  for (j in 1:nj) {
    for (i in 1:ni) {
      dist <- sqrt( (x.grd[i]-x)^2 + (y.grd[j]-y)^2 )
      i.near <- order(dist)
      i.near <- i.near[1:n.nearest]
      i.ok <- dist[i.near] <= dx*max.dist
      if (sum(i.ok) > 0) {
        i.near <- i.near[i.ok]
        if (min(dist[i.near]) > 0) {
          map.xy[i,j] <-  sum(map[i.near]/dist[i.near])/sum(1/dist[i.near])
        } else map.xy[i,j] <- map$map[i.near[1]]
      }
    }
  }

  z.levs <- seq(-max(abs(as.vector(map.xy)),na.rm=T),
                 max(abs(as.vector(map.xy)),na.rm=T),length=41)
  my.col <- rgb(c(seq(0,1,length=20),rep(1,21)),
                c(abs(sin((0:40)*pi/40))),
                c(c(rep(1,21),seq(1,0,length=20))))
  if (!add) {
    filled.contour(x.grd,y.grd,map.xy,
                   col = my.col,levels=z.levs,
                   main=paste(attributes(map.obj)$"long_name",
                              attributes(map.obj)$"descr"),
                   sub=date,xlab="",ylab="")
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

  points(x,y,pch=".")
# Circumpherence of the Earth:
  
  x <- seq(-1,1,length=100)
  y <- sqrt(1 - x^2)
  lines(x,y,type="l")
  lines(x,-y,type="l")
  
  data(addland)
#  keep <- (lat.cont-lat.0 >= -90) & (lat.cont-lat.0 <= 90) &
#          (lon.cont-lon.0 >= -90) & (lon.cont-lon.0 <= 90)
  keep <- (lat.cont >= 0)
  lat.cont[!keep] <- NA
  lon.cont[!keep] <- NA
  r.cont <- sin( pi*(lat.cont-lat.0)/180 )
  x.cont <- r.cont*sin(-pi*(lon.cont-lon.0)/180)
  y.cont <- r.cont*cos(-pi*(lon.cont-lon.0)/180)
  lines(x.cont,y.cont,col="grey30")
  contour(x.grd,y.grd,map.xy,add=T,lwd=lwd,lty=lty,col=col)
}
