satellite <- function(map.obj,col="black",lwd=2,lty=1,add=FALSE,
                      las = 1,lon.0=NULL,lat.0=NULL,method="normal",
                      ni=100,nj=100, n.nearest=4,max.dist=3,landdata="addland2") {
  if (class(map.obj)!="map")  stop("Need a map object (mapField)")

  if (!is.null(map$map.1)) n.maps <- map$n.maps else n.maps <- 1
  
  map.xy <- matrix(rep(NA,ni*nj*n.maps),ni,nj*n.maps)
  dim(map.xy) <- c(ni,nj,n.maps)
  for (i.map in 1:n.maps) {
    if (!is.null(map$n.maps)) {
      expression <- paste("map$map.",i.map,",[i.near[1]]",sep="")
      lon <- eval(parse(text=paste("map.obj$lon.",i.map,sep="")))
      lat <- eval(parse(text=paste("map.obj$lat.",i.map,sep="")))
    } else {
      lon <- map.obj$lon; lat <- map.obj$lat
      expression <- "map$map[i.near[1]]"
    }
 
    if (is.null(lon.0)) lon.0 <- mean(lon)
    if (is.null(lat.0)) lat.0 <- mean(lat)
  
    nx <- length(lon)
    ny <- length(lat)
    np <- nx*ny
    theta <- pi*lon.0/180
    phi <- pi*lat.0/180
    phi.min <- pi*min(lat[is.finite(lat)])/180
    r0 <- c(cos(phi)*cos(theta),
            sin(phi),
            cos(phi)*sin(theta))
    lats <- rep(lat,nx)
    lons <- sort(rep(lon,ny))

#    print(c(lon.0,lat.0))
#    print(range(lons))
#    print(range(lats))
#    print(r0)
  
    x <- rep(NA,np)
    y <- rep(NA,np)
    if (method=="polarstereo") {
      r <- sin( pi*(90-lats)/180 )
      x <- r*sin(pi*lons/180)
      y <- -r*cos(pi*lons/180)
    } else if (method=="distance") {
    for (i in 1:np) {
        theta <- pi*lons[i]/180
        phi <- pi*lats[i]/180
        r <- c(cos(phi)*cos(theta),
               sin(phi),
               cos(phi)*sin(theta))
        theta <- pi*lon.0/180
        r1 <- c(cos(phi)*cos(theta),
                sin(phi),
                cos(phi)*sin(theta))
        if (is.finite(r0*r1)) y[i] <- acos(sum(r0*r1)) else y[i] <- NA
        if (is.finite(r*r1)) x[i] <- acos(sum(r*r1)) else x[i] <- NA
      }
      y[y>0.5*pi] <- NA; x[x>0.5*pi] <- NA
      x <- sin(x); y <- sin(y)
      x[lons<lon.0] <- x[lons<lon.0]*-1
      y[lats<lat.0] <- y[lats<lat.0]*-1
      x[!is.finite(x)] <- 0
      y[!is.finite(y)] <- 0
    } else {
        for (i in 1:np) {
        theta <- pi*lons[i]/180
        phi <- pi*lats[i]/180
        r <- c(cos(phi)*cos(theta),
               sin(phi),
               cos(phi)*sin(theta))
        theta <- pi*lon.0/180
        phi <- phi.min
        r1 <- c(cos(phi)*cos(theta),
                sin(phi),
                cos(phi)*sin(theta))
        a <- r - r0
        b <- r1 - r0
        if (is.finite(r*r0)) newphi <- acos( sum(r*r0) ) else newphi<-NA 
        newtheta <- acos(sum(a*b) / (sqrt(sum(a*a)) * sqrt(sum(b*b))) )
        if (lons[i]<lon.0) newtheta <- -newtheta
        d <- sin(newphi)
        x[i] <- d*sin(newtheta)
        y[i] <- -d*cos(newtheta)
      }
    }

  
    map <- t(map.obj$map)
    dim(map) <- c(np,1)

    x.grd <- seq(-1,1,length=ni)
    y.grd <- seq(-1,1,length=nj)
    dx <- sqrt(mean(diff(x.grd))^2 + mean(diff(y.grd))^2)

#  print("gridding...")
#  x11(); plot(x,y)

    for (j in 1:nj) {
      for (i in 1:ni) {
        dist <- sqrt( (x.grd[i]-x)^2 + (y.grd[j]-y)^2 )
        i.near <- order(dist)
        i.near <- i.near[1:n.nearest]
        i.ok <- (dist[i.near] <= dx*max.dist) & (is.finite(map[i.near]))
        if ((sum(i.ok) > 0) & (x.grd[i]^2 + y.grd[j]^2 <= 1)) {
          i.near <- i.near[i.ok] 
          if (min(dist[i.near]) > 0) {
            map.xy[i,j,i.map] <-  sum(map[i.near]/dist[i.near])/sum(1/dist[i.near])
          } else map.xy[i,j,i.map] <- eval(parse(text=expression))
        } else map.xy[i,j,i.map] <- NA
      }
    }
  }

  z.levs <- seq(-max(abs(as.vector(map.xy)),na.rm=T),
                 max(abs(as.vector(map.xy)),na.rm=T),length=41)
  my.col <- rgb(c(seq(0,1,length=20),rep(1,21)),
                c(abs(sin((0:40)*pi/40))),
                c(c(rep(1,21),seq(1,0,length=20))))

  if (!add) {
    par(col.lab="white")
    filled.contour(x.grd,y.grd,map.xy[,,1],
                   col = my.col,levels=z.levs,
                   main=paste(attributes(map.obj)$"long_name",
                              attributes(map.obj)$"descr"),
                   sub=map.obj$date,xlab="",ylab="")
    par(col.lab="black")
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

#  points(x,y,pch=".")

# Circumpherence of the Earth:
#  

  if (!add) {
    eval(parse(text=paste("data(",landdata,")",sep="")))
    lon.cont[lon.cont > 180] <- lon.cont[lon.cont > 180] - 360

    iwest <- lon.cont < lon.0
    isouth <- lat.cont < lat.0
    np.cont <- length(lon.cont)
    y.cont <- rep(NA,np.cont)
    x.cont <- y.cont
    if (method=="polarstereo") {
      x <- seq(-1,1,length=100)
      y <- sqrt(1 - x^2)
      lines(x,y,type="l")
      lines(x,-y,type="l")

      r <- sin( pi*(90-lat.cont)/180 )
      x.cont <- r*sin(pi*(lon.cont)/180)
      y.cont <- -r*cos(pi*(lon.cont)/180)
      y.cont[lat.cont < 0] <- NA
    } else if (method=="distance") {
      for (i in 1:np.cont) {
        theta <- pi*lon.cont[i]/180
        phi <- pi*lat.cont[i]/180
        r <- c(cos(phi)*cos(theta),
               sin(phi),
               cos(phi)*sin(theta))
        theta <- pi*lon.0/180
        r1 <- c(cos(phi)*cos(theta),
                sin(phi),
                cos(phi)*sin(theta))
        y.cont[i] <- acos(sum(r0*r1))
        x.cont[i] <- acos(sum(r*r1))
      }
      y.cont[abs(y.cont)>0.5*pi] <- NA
      x.cont[abs(x.cont)>0.5*pi] <- NA
      x.cont <- sin(x.cont); y.cont <- sin(y.cont)
      x.cont[iwest] <- -x.cont[iwest]
      y.cont[isouth] <- -y.cont[isouth]
      x.cont[!is.finite(x.cont)] <- NA
      y.cont[!is.finite(y.cont)] <- NA

  } else {
      x <- seq(-1,1,length=100)
      y <- sqrt(1 - x^2)
      lines(x,y,type="l")
      lines(x,-y,type="l")
      for (i in 1:np.cont) {
        if (is.finite(lon.cont[i])) {

          theta <- pi*lon.cont[i]/180
          phi <- pi*lat.cont[i]/180
          r <- c(cos(phi)*cos(theta),
                 sin(phi),
                 cos(phi)*sin(theta))
          theta <- pi*lon.0/180
          r1 <- c(cos(phi)*cos(theta),
                  sin(phi),
                  cos(phi)*sin(theta))
          a <- r - r0
          b <- r1 - r0

          if (is.finite(r*r0)) newphi <- acos( sum(r*r0) ) else newphi<-NA  
          newtheta <- acos(sum(a*b) / (sqrt(sum(a*a)) * sqrt(sum(b*b))) )

          if (!is.finite(newphi))  newphi<- NA
          if (!is.finite(newtheta))  newtheta<- NA
          if (newphi > 0.5*pi) newphi<- NA
          if ((lon.cont[i]<lon.0) & (phi<pi*lat.0/180)) newtheta <- -newtheta
          if ((lon.cont[i]>lon.0) & (phi>pi*lat.0/180)) newtheta <- pi - newtheta
          if ((lon.cont[i]<lon.0) & (phi>pi*lat.0/180)) newtheta <- newtheta - pi
          d <- sin(newphi)
          x.cont[i] <- d*sin(newtheta)
          y.cont[i] <- -d*cos(newtheta)
        }     
      }
    }
  }
  lines(x.cont,y.cont,col="grey30")
  for (i.map in 1:n.maps){
    contour(x.grd,y.grd,as.matrix(map.xy[,,i.map]),add=T,lwd=lwd,lty=lty,col=col)
  }
} 
