meanField <- function(x,lon.rng=NULL,lat.rng=NULL,t.rng=NULL) {

  dat <- x$dat
  lons <- x$lon
  lats <- x$lat
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
  nj <- length(lats)
  ni <- length(lons)
  map <- matrix(rep(NA,ni*nj),ni,nj)
  for (j in 1:nj) {
    for (i in 1:ni) {
      map[i,j] <- mean(dat[,j,i],na.rm=T)
    }
  }
  
  results <- list(lon=lons,lat=lats,map=map,
                  tim=NULL,date=t.rng)
  class(results) <- "map"
#  attr(results) <- attr(x)
  attr(results,"descr") <- "Mean values"
  invisible(results) 
}
