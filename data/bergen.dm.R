obs0 <- read.table("bergen_florida_dm.txt",
                   col.names=c("station.number","day","month","year","precip","t2m"))
obs <- list(dd=obs0$day,mm=obs0$month,yy=obs0$year,
            t2m=obs0$t2m,precip=obs0$precip,
            obs.name=c("Daily mean temperature","Daily precipitation"),
            unit=c("deg C","mm/day"),station=obs0$station.number[1],
            lat=60.38,lon=5.33,alt=12,ele=c("tam","rr"),
            x.0E65N=NULL,y.0E65N=NULL,found=T,
            location="Bergen-Florida",wmo.no=NULL,
            start=1949,yy0=1980, country="Norway",ref="www.met.no")
class(obs) <- c("station","daily.station.record")
bergen.dm <- obs
rm(obs0)

# Stnr Navn I drift fra I drift til Snr Utm_e Utm_n Utm_sone Hoh Kommune Fylke # 50540 BERGEN - FLORIDA 1949   317 -31614 6733223 33 12 BERGEN HORDALAND
