# This function converts a time series to the same format as getnacd
# and getnordklim - compatible with ds.R
#
#
# R.E. Benestad 04.06.2002

station.obj <- function(x,yy,obs.name,unit,ele=NULL,mm=NULL,
                        station=NULL,lat=NULL,lon=NULL,alt=NULL,
                        location="unspecified",wmo.no=NULL,
                        start=NULL,yy0=NULL,country=NULL,ref=NULL) {

  if ((!is.null(lat)) & (!is.null(lon))) {
    xy <- COn0E65N(lon,lat)
 } else xy <- list(x=NULL,y=NULL)

  x[x <= -999] <-  NA
  yrs <- table(yy)
  ny <- length(row.names(yrs))
  yrs <- as.numeric(row.names(yrs))
  if (is.null(mm)) mm <- rep(1:12,ny)
#  print(table(yy))
#  print(table(mm))
#  print(c(mm[1],mm[length(mm)],ny))  
#  print(yrs)

#x11()
#plot(c(1,12),c(260,295),type="n")
#grid()
#col<-c("black","red","blue","green","darkred","darkblue","darkgreen","grey30")
#lwd <- c(2,rep(1,7),2)

#if (!is.null(mm)) print(paste("First month:",mm[1]," - last month:",mm[length(mm)]))  
  if (is.vector(x)) {
#    print("Vector")
    if (length(x)==length(yy)) {
      x.2D <- matrix(rep(NA,ny*12),ny,12)
      for (i in 1:ny) {
        x.i <- x[yy==yrs[i]]
        m.i <- mm[yy==yrs[i]]
 #       print(rbind(c(i,NA,m.i),c(yrs[i],NA,x.i)))
        x.2D[i,m.i] <- x.i
 #       lines(as.vector(x.2D[i,]),col=col[i],lwd=lwd[i])
 #       points(as.vector(x.2D[i,]),col=col[i],pch=20)
      }
    } else {
      print(paste("Different lengths: x->",length(x),
                  "    yy->",length(yy)))
      stop("Error: time stamp doesn't agree with data")
    }
    x <- x.2D
    yy <- yrs
    rm(x.2D,yrs); gc(reset=TRUE)
  }

#  print("Make list")
  station.obj<-list(val=x,yy=yy,station=station,
                    lat=lat,lon=lon,alt=alt,
                    x.0E65N=xy$x,y.0E65N=xy$y,
                    location=location, wmo.no=wmo.no,
                    start=start,yy0=NULL,ele=ele,
                    obs.name=obs.name, unit=unit,country=country,
                    found=TRUE,
                    ref=ref)
  class(station.obj) <- c("station","monthly.station.record")
  invisible(station.obj )
}
