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
#  print(table(yy))
#  print(table(mm))
#  print(c(mm[1],mm[length(mm)]))  
#  print(yrs)
  
  if (is.vector(x)) {
#    print("Vector")
    if (length(x)==length(yy)) {
      x.2D <- matrix(rep(NA,ny*12),ny,12)
      for (i in 1:ny) {
        x.i <- x[yy==yrs[i]]
        m.i <- mm[yy==yrs[i]]
        print(c(m.i,x.i))
        x.2D[i,m.i] <- x.i
      }
    }
    x <- x.2D
    yy <- yrs
    rm(x.2D,yrs)
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
  station.obj 
}
