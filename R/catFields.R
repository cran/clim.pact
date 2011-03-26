catFields <- function(field.1,field.2=NULL,lat=NULL,lon=NULL,
                      plot.interp=FALSE,interval.1=NULL,
                      interval.2=NULL,mon=NULL,demean=TRUE,silent=FALSE,
                      fastregrid=FALSE,neofs=20) {
  library(akima)
  l.one=FALSE
  l.newgrid <- FALSE
  if (is.null(field.2)) {
    l.one <- TRUE
    field.2 <- field.1
  }

  if ( (is.null(lat)) & (is.null(lon)) &
     (  (min(field.2$lon) > min(field.1$lon)) |
        (max(field.2$lon) < min(field.1$lon)) |
        (min(field.2$lat) > min(field.1$lat)) |
        (min(field.2$lat) < min(field.1$lat))  ) ) {
    if (!silent) print("Second field covers smaller domain:")
    i.lat <- (field.1$lat >= min(field.2$lat)) &
             (field.1$lat <= max(field.2$lat))
    i.lon <- (field.1$lon >= min(field.2$lon)) &
             (field.1$lon <= max(field.2$lon))
    field.1$lat <- field.1$lat[i.lat]
    field.1$lon <- field.1$lon[i.lon]
    field.1$dat <- field.1$dat[,i.lat,i.lon]
    field.1$id.lon <- field.1$id.lon[i.lon]
    field.1$id.lat <- field.1$id.lat[i.lat]
    field.1$id.x <- field.1$id.x[i.lat,i.lon]
    if (!silent) {print(range(field.1$lon)); print(range(field.1$lat))}
  }
  
  if (length(class(field.1))==2) {
    if (class(field.1)[1]=="mix.fields") {
      stop("Call mix.fields after catFields")
    }
  }
  if (class(field.1)[1] != class(field.2)[1]) {
    print(class(field.1)[1])
    print(class(field.2)[1])
    stop("The objects must have the same class")
  }

  tim.unit1 <- attr(field.1$tim,"unit")
  tim.torg1 <- attr(field.1$tim,"time_origin")
  tim.unit2 <- attr(field.2$tim,"unit")
  tim.torg2 <- attr(field.2$tim,"time_origin")
  if (is.null(tim.unit1)) tim.unit1<- "month"
  if (is.null(tim.unit2)) tim.unit2<- "month"
  if (lower.case(substr(tim.unit1,1,3)) != lower.case(substr(tim.unit2,1,3))) {
    print(c(tim.unit1,tim.unit2))
    stop('The time units must match')
  }
  
  if (!is.null(interval.1)) {
    if (!silent) print(interval.1)
    i1 <- ( (field.1$yy>=interval.1[1]) & (field.1$yy<=interval.1[2]))
    field.1$dat <- field.1$dat[i1,,]
    field.1$tim <- field.1$tim[i1]
    field.1$id.t <- field.1$id.t[i1]
    field.1$yy <- field.1$yy[i1]
    field.1$mm <- field.1$mm[i1]
    field.1$dd <- field.1$dd[i1]
  } else i1 <- is.finite(field.1$yy)
  
  if (!is.null(interval.2)) {
    if (!silent) print(interval.2)
    i2 <- ( (field.2$yy>=interval.2[1]) & (field.2$yy<=interval.2[2]))
    field.2$dat <- field.2$dat[i2,,]
    field.2$tim <- field.2$tim[i2]
    field.2$id.t <- field.2$id.t[i2]
    field.2$yy <- field.2$yy[i2]
    field.2$mm <- field.2$mm[i2]
    field.2$dd <- field.2$dd[i2]
  } else i2 <- is.finite(field.2$yy)
  
  if (!is.null(mon)) {
    cmon<-c('Jan','Feb','Mar','Apr','May','Jun',
            'Jul','Aug','Sep','Oct','Nov','Dec')
    if (!silent) print(paste("Extract",cmon[mon]))
    i1 <- is.element(field.1$mm,mon)
    i2 <- is.element(field.2$mm,mon)
    field.1$dat <- field.1$dat[i1,,]
    field.1$tim <- field.1$tim[i1]
    field.1$id.t <- field.1$id.t[i1]
    field.1$yy <- field.1$yy[i1]
    field.1$mm <- field.1$mm[i1]
    field.1$dd <- field.1$dd[i1]
    field.2$dat <- field.2$dat[i2,,]
    field.2$tim <- field.2$tim[i2]
    field.2$id.t <- field.2$id.t[i2]
    field.2$yy <- field.2$yy[i2]
    field.2$mm <- field.2$mm[i2]
    field.2$dd <- field.2$dd[i2]
  }
  nt.1 <- length(field.1$tim)
  nx.1 <- length(field.1$lon)
  ny.1 <- length(field.1$lat)
  nt.2 <- length(field.2$tim)
  nx.2 <- length(field.2$lon)
  ny.2 <- length(field.2$lat)
  
  #print(paste("Field1: ",nt.1,nx.1,ny.1,"   Field2: ",nt.2,nx.2,ny.2))
  
  if (xor(min(field.1$lon)<0,min(field.2$lon)<0)) {
    if (min(field.2$lon)<0) {
      field.1$lon[field.1$lon > 180] <- field.1$lon[field.1$lon > 180]-360
      x.srt <- order(field.1$lon)
      field.1$lon <- field.1$lon[x.srt]
      field.1$dat <- field.1$dat[,,x.srt]
    } else {
      field.2$lon[field.2$lon > 180] <- field.2$lon[field.2$lon > 180]-360
      x.srt <- order(field.2$lon)
      field.2$lon <- field.2$lon[x.srt]
      field.2$dat <- field.2$dat[,,x.srt]
    }
  }
  field.1$dat[!is.finite(field.1$dat)] <- 0
  field.2$dat[!is.finite(field.2$dat)] <- 0

  if ( (length(dim(field.1$dat)) !=3) |  length(dim(field.2$dat))!=3 ) {
      print("---------------A potential problem is detected in catFields !------------------")
      print(paste("dim(field.1$dat): sum(i1)=", sum(i1),"length(i1)=",length(i1)))
      print(dim(field.1$dat)); print(c(nt.1,ny.1,nx.1))
      print(paste("dim(field.2$dat): sum(i2)=",sum(i2),"length(i2)=",length(i2)))
      print(dim(field.2$dat)); print(c(nt.2,ny.2,nx.2))
      dim(field.1$dat) <- c(nt.1,ny.1,nx.1)
      dim(field.2$dat) <- c(nt.2,ny.2,nx.2)
  }

  if ((demean) & !l.one) {
    if (!silent) print("Subtracting mean values...")
#    dim(field.1$dat) <- c(nt.1,ny.1*nx.1)
#    field.1$dat <- field.1$dat - colMeans(field.1$dat)
#    dim(field.1$dat) <- c(nt.1,ny.1,nx.1)
#    dim(field.2$dat) <- c(nt.2,ny.2*nx.2)
#    field.2$dat <- field.2$dat - colMeans(field.2$dat)
#    dim(field.2$dat) <- c(nt.2,ny.2,nx.2)
### Slow...

    for (j in 1:ny.1) {
        for (i in 1:nx.1) {
          ave1 <- mean(field.1$dat[,j,i],na.rm=TRUE)
          if (!is.finite(ave1)) ave1 <- 0
          field.1$dat[,j,i] <- field.1$dat[,j,i]- ave1
        }
    }
         
    for (j in 1:ny.2) {
      for (i in 1:nx.2) {
        ave2 <- mean(field.2$dat[,j,i],na.rm=TRUE)
        if (!is.finite(ave2)) ave2 <- 0
        field.2$dat[,j,i] <- field.2$dat[,j,i]- ave2
      }
    }
    if (!silent) print("---")
  }
  if (!is.null(lat) & !is.null(lon)) {
    interpolate <- !((length(lat)==2) & (length(lon)==2))
    if (length(lat)==2) {
      iy <- field.1$lat >= min(lat) & field.1$lat <= max(lat)
      lat <- field.1$lat[iy]
    }
    if (length(lon)==2) {
      ix <- field.1$lon >= min(lon) & field.1$lon <= max(lon)
      lon <- field.1$lon[ix]
    }
    ny.1 <- length(lat); nx.1 <- length(lon)
    field.1$id.x <-  matrix(rep(field.1$v.nam,ny.1*nx.1),ny.1,nx.1)
    field.1$id.lon <- rep(field.1$v.nam,nx.1); field.1$id.lat <- rep(field.1$v.nam,ny.1)

    if (!interpolate) {
      #print(dim(field.1$dat))
      dat.1 <- field.1$dat[,iy,ix]
      #print(dim(dat.1))
    } else {
      l.newgrid <- TRUE
      if (!fastregrid) {
        if (!silent) print("interpolate 1st field - please be patient")
        lat.x<-rep(field.1$lat,length(field.1$lon))
        lon.x<-sort(rep(field.1$lon,length(field.1$lat)))
        dat.1<-matrix(nrow=nt.1,ncol=ny.1*nx.1)
        dim(dat.1)<-c(nt.1,ny.1,nx.1)
        for (it in 1:nt.1) {
            Z.in<-as.matrix(field.1$dat[it,,])
            Z.out<-interp(lat.x,lon.x,Z.in,lat,lon)
            dat.1[it,,]<-as.matrix(Z.out$z)
            if (plot.interp) {
               contour(field.1$lon,field.1$lat,t(round(Z.in,2)),col="blue",lwd=2,
                       main=paste("Field 1: ",it,"/",nt.1),
                       sub=paste(field.1$mm[it],"-",field.1$yy[it]))
              contour(lon,lat,t(round(Z.out$z,2)),add=TRUE,col="red",lty=2)
              addland()
              grid()
            }
          }
        } else {   # REB 02.12.2010 - new lines to speed up interpolation
          if (!silent) print("First field - quick interpolation :-)")
          regridded <- fastRegrid(field.1,lat=lat,lon=lon,neofs=neofs)
          dat.1 <- regridded$dat
          lat <- regridded$lat
          lon <- regridded$lon
          nx.1 <- length(regridded$lon)
          ny.1 <- length(regridded$lat)
      }
    }
  } else {
    if (!silent) print("Use the grid of first field - no interpolation :-)")
    lat <- field.1$lat
    lon <- field.1$lon
    dat.1 <- field.1$dat
    nx.1 <- length(lon)
    ny.1 <- length(lat)
  }
  
  lat.x<-rep(field.2$lat,length(field.2$lon))
  lon.x<-sort(rep(field.2$lon,length(field.2$lat)))
  dat.2<-matrix(nrow=nt.2,ncol=ny.1*nx.1)
  dim(dat.2)<-c(nt.2,ny.1,nx.1)
  l.different <- TRUE
  if ( (ny.1==ny.2) & (nx.1==nx.2) ) {
#    if ( (sum(field.1$lat==field.2$lat)==ny.1) &
#         (sum(field.1$lon==field.2$lon)==nx.1) ) l.different <- FALSE
# REB 29.03.2005: problem with above lines: 'Error in field.1$lat == field.2$lat : non-conformable arrays'
    if ( (sum(is.element(field.1$lat,field.2$lat))==ny.1) &
         (sum(is.element(field.1$lon,field.2$lon))==nx.1) ) l.different <- FALSE
  }

#  print(c(l.one,l.different,l.newgrid))

  if (nt.2==0) {
    print("catFields: nt.2=0!")
    print(dim(field.2$dat))
    return(NULL)
  }
  if (!l.one & (l.different | l.newgrid)) {
    if (!fastregrid) {
      if (!silent) print("Interpolate 2nd field - please be patient")
      for (it in 1:nt.2) {
        Z.in<-as.matrix(field.2$dat[it,,])
        Z.out<-interp(lat.x,lon.x,Z.in,lat,lon)
        dat.2[it,,]<-as.matrix(Z.out$z)
        if (plot.interp) {
          contour(field.2$lon,field.2$lat,t(round(Z.in,2)),col="blue",lwd=2,
                  main=paste("Field 1: ",it,"/",nt.2),
                  sub=paste(field.2$mm[it],"-",field.2$yy[it]))
          contour(lon,lat,t(round(Z.out$z,2)),add=TRUE,col="red")
          addland()
          grid()
        }
      }
    } else {   # REB 02.12.2010 - new lines to speed up interpolation
      if (!silent) print("Interpolate 2nd field - using fastRegrid")
      regridded <- fastRegrid(field.2,lat=lat,lon=lon,mon=mon,neofs=neofs)
      dat.2 <- regridded$dat
      nt.2 <- length(regridded$tim)
      #print(length(regridded$tim));print(dim(dat.2));print(c(nt.2,ny.1,nx.1))
    }
  } else {
    #print("Identical spatial grids :-)")
    dat.2 <- field.2$dat
  }
  #print(dim(dat.1));  print(dim(dat.2)); print(c(nt.1,ny.1*nx.1));
  #print(c(nt.2,ny.2*nx.2)); print(table(field.1$mm)); print(table(field.2$mm))
  
  dim(dat.1)<-c(nt.1,ny.1*nx.1)
  if (!l.one) {
    dim(dat.2)<-c(nt.2,ny.1*nx.1)
    dat<-rbind(dat.1,dat.2)
    dim(dat)<-c(nt.1+nt.2,ny.1,nx.1)
    tim <- c(field.1$tim,field.2$tim)
    yy <- c(field.1$yy,field.2$yy)
    mm <- c(field.1$mm,field.2$mm)
    dd <- c(field.1$dd,field.2$dd)
    id.t <- c(field.1$id.t,field.2$id.t)
  } else {
    dat<-dat.1
    dim(dat)<-c(nt.1,ny.1,nx.1)
    tim <- field.1$tim
    yy <- field.1$yy
    mm <- field.1$mm
    dd <- field.1$dd
    id.t <- field.1$id.t
  }
  id.x <- matrix(rep(field.1$id.x[1],ny.1*nx.1),ny.1,nx.1)
  attr(tim,"unit") <- tim.unit1
  attr(tim,"time_origin") <- tim.torg1
  if (field.1$v.name==field.2$v.name) var.name <- field.1$v.name else
  var.name <- paste(field.1$v.name,"&",field.2$v.name,sep="")
  result  <- list(dat=dat,lon=lon,lat=lat,tim=tim,v.name=var.name,
                  id.t=id.t,id.x=id.x,yy=yy,mm=mm,dd=dd,n.fld=field.1$n.fld,
                  id.lon=field.1$id.lon,id.lat=field.1$id.lat,attributes=field.1$attributes,
                  filename=paste(field.1$filename,"+",field.2$filename))
  class(result) <- c(class(field.1),"cat.fields")
  invisible(result)
}

fastRegrid <- function(field.1,lat=NULL,lon=NULL,mon=NULL,neofs=20,
                       silent=TRUE,plot.interp=FALSE) {
  fieldIDs <- rownames(table(field.1$id.x))
  if (length(fieldIDs)>1) {
    print("Can only handle uni-fields")
    print(table(field.1$id.x))
    return()
  }
  if (is.null(lat)) lat <- field.1$lat
  if (is.null(lon)) lon <- field.1$lon
  nx <- length(lon); ny <- length(lat)
  mean.field <- field.1$dat
  d <- dim(mean.field)
  dim(mean.field) <- c(d[1],d[2]*d[3])
  mean.field <- colMeans(mean.field,na.rm=TRUE)
  
  eof <- EOF(field.1,lon=range(lon),lat=range(lat),mon=mon,
             neofs=neofs,plot=plot.interp,silent=silent)
  lat.x<-rep(eof$lat,length(eof$lon))
  lon.x<-sort(rep(eof$lon,length(eof$lat)))
  regridd.eof <- rep(NA,nx*ny*neofs)
  dim(regridd.eof) <- c(neofs,nx*ny)
  for (it in 1:neofs) {
        Z.in<-as.matrix(eof$EOF[it,])
        dim(Z.in) <- c(eof$size[2:3])
        good <- is.finite(Z.in)
        Z.out<-as.matrix(interp(lat.x[good],lon.x[good],Z.in[good],lat,lon)$z)
        dim(Z.out) <- c(ny*nx)
        regridd.eof[it,] <- Z.out
  }
  Z.in<-as.matrix(eof$clim)
  good <- is.finite(Z.in)
  Z.mean<-as.matrix(interp(lat.x[good],lon.x[good],Z.in[good],lat,lon)$z)
  eof$clim <- Z.mean
  eof$EOF <- regridd.eof
  eof$size[2] <- ny
  eof$size[3] <- nx
  eof$lon <- lon
  eof$lat <- lat
  eof$id.lon <- rep(eof$id.lon[1],nx)
  eof$id.lat <- rep(eof$id.lat[1],ny)
  eof$id.x <- c(eof$id.lat,eof$id.lon)
  field.2 <- EOF2field(eof)
  field.2$history <- "fastRegrid"
  invisible(field.2)
}

testfastRegrid <- function() {
  data(oslo.t2m)
  data(DNMI.t2m)
  data(DNMI.slp)
  slp <- fastRegrid(DNMI.slp,lat=seq(30,80,by=2.5),lon=seq(-90,60,by=2.5))
  mapField(DNMI.slp)
  mapField(slp)
  slp1 <- plotField(slp,lon=0,lat=60,what="abs")
  slp2 <- plotField(DNMI.slp,lon=0,lat=60,what="abs",add=TRUE,col="red",lty=2)

  newFig()
  eof1 <- EOF(slp,mon=1); newFig()
  eof2 <- EOF(DNMI.slp,mon=1)

  plotEOF(eof1)
  plotEOF(eof2)
  
  corField(DNMI.slp,oslo.t2m,mon=1); newFig()
  corField(slp,oslo.t2m,mon=1)

  t2m <- fastRegrid(DNMI.t2m,lat=seq(30,80,by=1),lon=seq(-90,60,by=1))

  newFig()
  plotStation(oslo.t2m,l.anom=FALSE,what="t",mon=7,pch=19,
              col="grey",ylim=c(5,25))
  t2m1 <- plotField(t2m,lon=oslo.t2m$lon,lat=oslo.t2m$lat,
                    what="abs",add=TRUE,col="blue",mon=7,lwd=2)
  t2m2 <- plotField(DNMI.t2m,lon=oslo.t2m$lon,lat=oslo.t2m$lat,
                    what="abs",add=TRUE,col="red",lty=2,mon=7)

  eof1b <- EOF(t2m,mon=1); newFig()
  eof2b <- EOF(DNMI.t2m,mon=1)
  
  ds1 <- DS(oslo.t2m,eof1b)
  ds2 <- DS(oslo.t2m,eof2b)
  
}
