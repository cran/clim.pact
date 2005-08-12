obs0 <- read.table("oslo_blindern_dm.txt",
                   col.names=c("station.number","day","month","year","t2m","precip"))
oslo.dm <- list(day=obs0$day,month=obs0$month,year=obs0$year,
            t2m=obs0$t2m,precip=obs0$precip,
            obs.name=c("Daily mean temperature","Daily precipitation"),
            unit=c("deg C","mm/day"),station=obs0$station.number[1],
            lat=67.33223,lon=-3.1614,alt=12,ele=c("tam","rr"),
            x.0E65N=NULL,y.0E65N=NULL,found=TRUE,
            location="Oslo-Blindern",wmo.no=NULL,
            start=1949,yy0=1980, country="Norway",ref="www.met.no")
class(oslo.dm) <- "daily.station.record"
rm(obs0)
# Stnr Navn I drift fra I drift til Snr Utm_e Utm_n Utm_sone Hoh Kommune Fylke # 50540 BERGEN - FLORIDA 1949   317 -31614 6733223 33 12 BERGEN HORDALAND
