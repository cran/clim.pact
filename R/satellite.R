
stereogr<- function(map.obj,NH=TRUE,lat.0=0,inv.col=FALSE,levels=NULL,sym=TRUE,dr=0.01,main=NULL) {
  old.par <- par()
  nc1 <- 20; nc2 <- 21
  if (sym) {
    z.levs <- seq(-max(abs(as.vector(map.obj$map)),na.rm=TRUE),
                   max(abs(as.vector(map.obj$map)),na.rm=TRUE),length=41)
  } else {
    z.levs <- seq(min(as.vector(map.obj$map),na.rm=TRUE),
                  max(as.vector(map.obj$map),na.rm=TRUE),length=41)
  }
  nc1 <- 20; nc2 <- 21
  if (!is.null(levels)) {
    z.levs <- levels
    nc1 <- floor(length(levels)/2)
    nc2 <- ceiling(length(levels)/2)
  }
  
  my.col <- rgb(c(seq(0,1,length=nc1),rep(1,nc2)),
                c(abs(sin((0:(length(z.levs)-1))*pi/(length(z.levs)-1)))),
                c(c(rep(1,nc2),seq(1,0,length=nc1))))
    if (inv.col) my.col <- reverse(my.col)
  par(xaxt="n",yaxt="n")
  nx <- length(map.obj$lon); ny <- length(map.obj$lat)
  lonx <- rep(map.obj$lon,ny); latx <- sort(rep(map.obj$lat,nx))
  if (!NH) latx <- -latx
  r <- sin( pi*(90-latx)/180 )
  x <- r*sin(pi*lonx/180)
  y <- -r*cos(pi*lonx/180)
  x.grd <- seq(-sin( pi*(90-lat.0)/180 ),sin( pi*(90-lat.0)/180 ),by=0.01);
  y.grd <- x.grd; Z <- c(map.obj$map)

  xylims <- max(r)*c(-1,1)
  nxy <- length(x.grd)
  good <- is.finite(Z) & (latx>lat.0)
  polar <- interp(x[good],y[good],Z[good],x.grd,y.grd,duplicate="mean")$z

  if (sum(!is.finite(map.obj$map))>0) {
    map <- polar+NA
    x.good <- x[good]; y.good <- y[good]
    X.grd <- rep(x.grd,nxy); Y.grd <- sort(rep(y.grd,nxy))
# test:     print(sum(good))
    for (i in 1:sum(good)) {
      r <- sqrt( (X.grd - x.good[i])^2 + (Y.grd - y.good[i])^2 )
      valid <- (r  < dr)
      map[valid] <- polar[valid]
    }
  } else map <- polar
# test:     print(summary(r)); print(length(r)); print(length(X.grd))
  if (is.null(main)) main=map.obj$description
  if (nchar(main)>40) par(cex.main=0.75)
  if (nchar(main)>60) par(cex.main=0.60)
  
  image(x.grd,y.grd,map,xlab="",ylab="",main=main,col = my.col,sub=map.obj$date,
        xlim=xylims,ylim=xylims)
  contour(x.grd,y.grd,map,add=TRUE)
# test:   points(x[good],y[good],pch="+",col="red")

  data(addland1,envir = environment())
  if (!NH) lat.cont <- -lat.cont
  ok <- is.finite(lon.cont) & is.finite(lat.cont)  & (lat.cont>lat.0)
  lat.cont <- abs(lat.cont)
  lon.cont[!ok] <-  -9999
  lat.cont[!ok] <-  -9999
  lon.cont[lon.cont > 180] <- lon.cont[lon.cont > 180] - 360
  lon.cont[!ok] <- NA; lat.cont[!ok] <- NA

  r <- sin( pi*(90-lat.cont)/180 )
  x.cont <- r*sin(pi*(lon.cont)/180)
  y.cont <- -r*cos(pi*(lon.cont)/180)
  lines(x.cont,y.cont,col="grey30")

  s <- seq(-2*pi,2*pi,length=360)
  lines(cos(s),sin(s))
  for (i in 1:8) {
    if (i*10 > lat.0) lines(sin( pi*(90-i*10)/180 )*cos(s),sin( pi*(90-i*10)/180 )*sin(s),col="grey")
    text(0,sin( pi*(90-i*10)/180 ),as.character(i*10),col="grey",cex=0.7)
  }
  stereog <- list(lon=x.grd,lat=y.grd,map=map,date=map.obj$date,tim=map.obj$tim,
                  description=map.obj$description)
  class(stereog) <- c("map","polar-stereographic")
  invisible(stereog)
#  par(old.par)
}

map2sphere <- function(z,X=seq(-180,180,by=1),Y=seq(-90,90,by=1),pal="rainbow",nlevs=21,breaks=NULL,quick=FALSE,add=FALSE) {
  print("map2sphere:")
  colour <- eval(parse(text=paste(pal,"(n=",nlevs,")",sep="")))
  if (class(z)=="map") {
    Z <- z$map
    X <- z$lon
    Y <- z$lat
    z <- Z
  } else if (is.list(z)) {
    components <- names(z)
    v1 <- eval(parse(text=paste("z$",components[1],sep="")))
    v2 <- eval(parse(text=paste("z$",components[2],sep="")))
    v3 <- eval(parse(text=paste("z$",components[3],sep="")))
    l <- c(length(v1),length(v2),length(v3))
    srt <- order(l,decreasing=TRUE)
    Z <- eval(parse(text=paste("z$",components[srt[1]],sep="")))
    d <- dim(Z)
    X <- eval(parse(text=paste("z$",components[l==d[1]],sep="")))
    Y <- eval(parse(text=paste("z$",components[l==d[2]],sep="")))
    X[X > 180] <- X[X > 180] - 360
    srtx <- order(X)
    X <- X[srtx]
    z <- Z[srtx,]
  }
  if (is.null(breaks)) breaks <- seq(min(z,na.rm=TRUE),max(z,na.rm=TRUE),length=nlevs+1)
  
  # print(X);print(Y)
  nx <- length(z[,1]); ny <-  length(z[1,])
  #print(c(nx,ny))
  nX <- length(X); nY <- length(Y)
  dx <- 0.5*diff(X)[1]; dy <- 0.5*diff(Y)[1]
  if ( (nx != length(X)) | (ny != length(Y)) ) {
    x <- rep((1:nx),ny); y <- sort(rep((1:ny),nx))
    Z <- interp(x,y,z,X,Y)$z
    z <- Z
  }
  x11()
  if (!add) plot(cos(seq(0,2*pi,length=360)),sin(seq(0,2*pi,length=360)),type="l",xlab="",ylab="")
  phi <- c(-80,-70,-60,-50,-40,40,50,60,70,80,90)
  steps <- c(6,5,4,3,2,1,2,3,4,5,6)
  j1 <- 1
  for (jj in 1:length(phi)) {
    j2 <- max( (1:ny)[Y <= phi[jj]] )
    ystep <- min(c(steps[jj],j2-j1))
    if (!quick)  ystep <- 1 
    sy1 <- sin(2*pi*(Y-dy*ystep)/360)
    sy2 <- sin(2*pi*(Y+dy*ystep)/360)
    cy <- cos(2*pi*Y/360)
    #print(c(j1,j2,ystep))
    for (j in seq(j1,j2,ystep)) {
      py1 <- sy1[j]
      py2 <- sy2[j]
      xstep <- steps[min( (1:length(steps))[phi >= Y[j]] )]
      if (!quick)  xstep <- 1
      x1 <- sin(2*pi*(X-dx*xstep)/360)*cy[j]
      x2 <- sin(2*pi*(X+dx*xstep)/360)*cy[j]
      #print(c(Y[j],xstep))
      for (i in seq(1,nX,by=xstep)) {
        if ( (X[i] >= -90) & (X[i] <= 90) & (is.finite(z[i,j])) ) {
          px1 <- x1[i] 
          px2 <- x2[i] 
          ilev <- round( nlevs*(z[i,j] - min(z,na.rm=TRUE))/diff(range(z,na.rm=TRUE)) )
          polygon(c(px2,rep(px1,2),rep(px2,2)),
                  c(rep(py1,2),rep(py2,2),py1),
                  col=colour[ilev], border=colour[ilev])
         }
       }
    }
    j1 <- j2+1
  }
  
}



satellite <- function(map,lon.0=0,lat.0=0,pal="rainbow",nlevs=21,breaks=NULL,quick=FALSE,add=FALSE) {
#  HERE
#  call map2sphere
#  Find map senter - flip if necessary
#  project onto sphere the visible part

  print("This function is not finished")
  map2sphere(z=map$map,X=map$lon,Y=map$lat,pal=pal,nlevs=nlevs,
             breaks=breaks,quick=quick,add=add)
}


satelliteOld <- function(map.obj,col="black",lwd=2,lty=1,add=FALSE,
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

  z.levs <- seq(-max(abs(as.vector(map.xy)),na.rm=TRUE),
                 max(abs(as.vector(map.xy)),na.rm=TRUE),length=41)
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
    eval(parse(text=paste("data(",landdata,",envir = environment())",sep="")))
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
    contour(x.grd,y.grd,as.matrix(map.xy[,,i.map]),add=TRUE,lwd=lwd,lty=lty,col=col)
  }
} 
