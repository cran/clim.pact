meanField <- function(x,lon.rng=NULL,lat.rng=NULL,t.rng=NULL,mon=NULL) {

  season<-cbind(c(12,1,2),c(3,4,5),c(6,7,8),c(9,10,11))
  months<-c("Jan","Feb","Mar","Apr","May","Jun",
          "Jul","Aug","Sep","Oct","Nov","Dec")
  season.c<-c("","DJF","MAM","JJA","SON")
  dat <- x$dat
  lons <- x$lon
  lats <- x$lat
  mm <- x$mm
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
    it <- (x$yy >= t.rng[1]) & (x$yy <= t.rng[2])
    dat <- dat[it,,]
  }
  if (!is.null(mon)) {
    if (class(x)[2]=="monthly.field.object") {
      i.mm <- is.element(mm,mon)
      dat <- dat[i.mm,,]
      mm <- mm[i.mm]
      month <- months[mon+1]
    } else if (class(x)[2]=="daily.field.object") {
      mon <- mod(mon-1,4)+1
      mon <- season[,mon]
      i.mm <- is.element(mm,mon)
      dat <- dat[i.mm,,]
      mm <- mm[i.mm]
      month <- season.c[mon+1]
    }
  }
  nj <- length(lats)
  ni <- length(lons)
  map <- matrix(rep(NA,ni*nj),ni,nj)
  for (j in 1:nj) {
    for (i in 1:ni) {
      map[i,j] <- mean(dat[,j,i],na.rm=T)
    }
  }
  
  results <- list(lon=lons,lat=lats,map=map,v.name=,x$v.name,
                  tim=month,date=t.rng)
  class(results) <- "map"
#  attr(results) <- attr(x)
  attr(results,"descr") <- "Mean values"
  invisible(results) 
}
