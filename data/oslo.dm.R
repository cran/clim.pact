obs0 <- read.table("oslo_blindern_dm.txt",
                   col.names=c("station.number","day","month",
                     "year","t2m","precip"))
oslo.dm <- list(dd=obs0$day,mm=obs0$month,yy=obs0$year,
            t2m=obs0$t2m,precip=obs0$precip,
            obs.name=c("Daily mean temperature","Daily precipitation"),
            unit=c("deg C","mm/day"),station=obs0$station.number[1],
            lat=59.95,lon=10.72,alt=94,ele=c("tam","rr"),
            x.0E65N=NULL,y.0E65N=NULL,found=T,
            location="Oslo-Blindern",wmo.no=NULL,
            start=1937,yy0=1980, country="Norway",ref="www.met.no")
class(oslo.dm) <- c("station","daily.station.record")
rm(obs0)

# Stnr Navn I drift fra I drift til Snr Utm_e Utm_n Utm_sone Hoh Kommune Fylke
# 18700 OSLO - BLINDERN feb 1937   492 261009 6652760 33 94 OSLO OSLO
