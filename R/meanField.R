
meanField <- function(x,lon.rng=NULL,lat.rng=NULL,t.rng=NULL,mon=NULL) {

  season<-cbind(c(12,1,2),c(3,4,5),c(6,7,8),c(9,10,11))
  months<-c("Jan","Feb","Mar","Apr","May","Jun",
            "Jul","Aug","Sep","Oct","Nov","Dec")
  season.c<-c("","DJF","MAM","JJA","SON")
  dat <- x$dat
  lons <- x$lon
  lats <- x$lat
  mm <- x$mm; yy <- x$yy; dd <- x$dd
  month <- NULL
  if (!is.null(lon.rng)) {
    ix <- (lons >= lon.rng[1]) & (lons <= lon.rng[2])
    dat <- dat[,,ix]
    lons <- lons[ix]
  }
  if (!is.null(lat.rng)) {
    iy <- (lats >= lat.rng[1]) & (lats <= lat.rng[2])
    dat <- dat[,iy,]
    lats <- lats[iy]
  }
  if (!is.null(t.rng)) {
    if (is.numeric(t.rng)) it <- (x$yy >= t.rng[1]) & (x$yy <= t.rng[2])
    if (is.character(t.rng)) {
      t.rng1 <- datestr2num(t.rng[1],vec=FALSE)
      t.rng2 <- datestr2num(t.rng[2],vec=FALSE)
      it <- (x$yy + x$mm/12 + x$mm/31 >= t.rng1) & 
            (x$yy + x$mm/12 + x$mm/31 <= t.rng2)
    }
    #print(c(sum(it),length(it),NA,length(dat[,1,1])))

    dat <- dat[it,,]
    yy <- yy[it]; mm <- mm[it]; dd <- dd[it]
  }

# Testing....
#  iy <- (1:length(lats))[lats >= 69.6538][1]
#  ix <- (1:length(lons))[lons >= 18.9283][1]
#  plot(yy+mm/12,dat[,iy,ix],ylim=c(-20,20),col="grey",type="b",pch=20)
#  print(range(t2m$dat[,iy,ix])); print(table(mm))
#  print(class(x)[2])

  if (!is.null(mon)) {
    if (class(x)[2]=="monthly.field.object") {
      i.mm <- is.element(mm,mon)
      dat <- dat[i.mm,,]
      yy <- yy[i.mm]
      mm <- mm[i.mm]
      dd <- dd[i.mm]
      month <- months[mon[1]]
      if (length(mon) > 1) 
        for (i in 2:length(mon)) month <- paste(month,"-",months[mon[i]],sep="")
    } else if (class(x)[2]=="daily.field.object") {
      mon <- mod(mon-1,4)+1
      mon <- season[,mon]
      i.mm <- is.element(mm,mon)
      dat <- dat[i.mm,,]
      yy <- yy[i.mm]
      mm <- mm[i.mm]
      dd <- dd[i.mm]
      month <- season.c[mon+1]
    }
  }
  nj <- length(lats)
  ni <- length(lons)
  nt <- length(mm)
  dim(dat) <- c(nt,nj*ni)
  map <- colMeans(dat)
  dim(dat) <- c(nt,nj,ni)
  dim(map) <- c(nj,ni)
  
# continued testing...
#  points(yy+mm/12,dat[,iy,ix])
#  lines(yy+mm/12,dat[,iy,ix])
#  x11()


#  map <- matrix(rep(NA,ni*nj),ni,nj)
#  for (j in 1:nj) {
#    for (i in 1:ni) {
#      map[i,j] <- mean(dat[,j,i],na.rm=T)
#    }
#  }
  
  results <- list(lon=lons,lat=lats,map=t(map),v.name=x$v.name,
                  tim=month,date=t.rng)
  class(results) <- "map"
#  attr(results) <- attr(x)
  attr(results,"long_name")<- paste("Mean",x$v.name)   
  attr(results,"descr") <- c("Mean values",month)
  invisible(results) 
}
