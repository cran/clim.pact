EOF2field <- function(eof,anomalies=FALSE) {
  x <-t(eof$EOF) %*% diag(eof$W) %*% t(eof$PC); x <- t(x)
  #print(dim(x)); print(eof$size)
  dim(x) <- eof$size
  if ((!anomalies) & !is.null(eof$clim)) { 
    clim <- t(matrix(eof$clim,eof$size[3],eof$size[2]))
    for (i in 1:eof$size[1]) x[i,,] <- x[i,,] + clim
  }
  retrieve.nc  <- list(dat=x,lon=eof$lon,lat=eof$lat,tim=eof$tim,lev=NULL,
                       v.name=eof$v.name,id.x=eof$id.x,id.t=eof$id.t,
                       yy=eof$yy,mm=eof$mm,dd=eof$dd,n.fld=eof$n.fld,
                       id.lon=eof$id.lon,id.lat=eof$id.lat,
                       attributes=eof$attributes)
  class(retrieve.nc) <- c(class(eof)[-(1)])
  invisible(retrieve.nc)
}
