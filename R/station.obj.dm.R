# This function converts a time series to the "daily.station.record"
# class - compatible with ds.dm.R
#
#
# R.E. Benestad 04.06.2002

station.obj.dm <- function(t2m,precip,dd,mm,yy,
                        obs.name=NULL,unit=NULL,ele=NULL,
                        station=NULL,lat=NULL,lon=NULL,alt=NULL,
                        location="unspecified",wmo.no=NULL,
                        start=NULL,yy0=NULL,country=NULL,ref=NULL) {

#  source("COn0E65N.R")
  xy<-COn0E65N(lon,lat)
#  if ((lat != NULL) & (lon != NULL)) xy<-COn0E65N(lon,lat) else
#     xy <- list(x=NULL,y=NULL)

  station.obj.dm<-list(t2m=t2m,precip=precip,
                    dd=dd,mm=mm,yy=yy,
                    obs.name=obs.name,unit=unit,ele=ele,
                    station=station,
                    lat=lat,lon=lon,alt=alt,
                    x.0E65N=xy$x,y.0E65N=xy$y,
                    location=location, wmo.no=wmo.no,
                    start=start,yy0=yy0,country=country,
                    found=TRUE,ref=ref)
  class(station.obj.dm) <- c("station","daily.station.record")
  station.obj.dm 
}
