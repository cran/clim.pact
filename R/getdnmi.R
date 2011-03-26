# This R routine reads the DNMI/met.no data. The code
# will not work with the original NACD files: a space
# must be inserted between the December value and the
# country tag, and the missing values must be changed
# from '-9999' to ' -999'.
#
# Arguments:
# 'location' determines the time series.
# 'ele.c' determines the element (default=T2m).
#
# R.E. Benestad

getdnmi <- function(location,ele.c='101',silent = FALSE,
                    direc="data/") {

if (!silent) print(paste("GETDNMI:",location,ele.c))
row<-switch(as.character(ele.c),
              't2m'='V5','rr'='V4','slp'='V9',
              't2'='V5','precip'='V4','temp'='V5',
              'snow depth'='V6','snowy days'='V7',
              'humidity'='V8','101'='V5','601'='V4',
              '112'='V10','122'='V11',
              'tax'='V10','tan'='V11')

descr <- switch(as.character(ele.c),
              't2m'='Temperature',
              'rr'='Precipitation',
              'slp'='Sea level pressure',
              't2'='Temperature',
              'precip'='Precipitation',
              'temp'='Temperature',
              'snow depth'='Snow depth',
              'snowy days'='Snowy days',
              'humidity'='Relative humidity',
              '101'='Temperature',
              '601'='Precipitation',
              '112'='Abs. max monthly Temp.',
              '122'='Abs. min monthly Temp.',
              'tax'='Abs. max monthly Temp.',
              'tan'='Abs. min monthly Temp.')
unit <- switch(as.character(ele.c),
              't2m'='deg C',
              'rr'='mm',
              'slp'='hPa',
              't2'='deg C',
              'precip'='mm',
              'temp'='deg C',
              'snow depth'='cm',
              'snowy days'='days',
              'humidity'='%',
              '101'='deg C',
              '601'='mm',
              '112'='deg C',
              '122'='deg C',
              'tax'='deg C',
              'tan'='deg C')

f.name<-paste(lower.case(location),".dat",sep="")
if (file.exists(paste(direc,f.name,sep=""))) {
  no.find <- FALSE 

  obs<-read.table(paste(direc,f.name,sep=""))
#  dnmi.meta <- read.table(paste(direc,"dnmi.meta",sep=""))
#  dnmi.meta <- read.table("data/dvh.station.list",header=TRUE,as.is=TRUE)
  data(dvh.station.list)
  station<-obs$V1[1]
#  alt<- dnmi.meta$alt[dnmi.meta$dnmi.no==station]
#  lon<- dnmi.meta$lon[dnmi.meta$dnmi.no==station]
#  lat<- dnmi.meta$lat[dnmi.meta$dnmi.no==station]
#  location <- dnmi.meta$location[dnmi.meta$dnmi.no==station]
  alt<- as.numeric(dnmi.meta$Hoh[dnmi.meta$Stnr==station][1])
  lon<- dnmi.meta$Lon[dnmi.meta$Stnr==station][1]
  lat<- dnmi.meta$Lat[dnmi.meta$Stnr==station][1]
  location <- dnmi.meta$Navn[dnmi.meta$Stnr==station][1]
  yy <- obs$V2
  mm <- obs$V3
  YY<-as.numeric(row.names(table(obs$V2)))
  ny <- length(YY)
  x <- rep(NA,ny*12)
  eval(parse(text=paste("dat <- obs$",row,sep="")))
  for (i in 1:(ny*12)) {
    i.yy <- yy[1] + floor(i/12)
    i.mm <- mod(i-1,12)+1
    ii <- (yy == i.yy) & (mm == i.mm)
    if (sum(ii)>0) {
       x[i] <- dat[ii]
    } else x[i] <- NA
  }
  x[x<= -999] <- NA
  dim(x) <- c(12,ny)
  x <- t(x)
  
  country<-"Norway"
  xy<-COn0E65N(lon,lat)
  yy<-row.names(table(obs$V2))
} else {
  x <- NULL
  station <- NA
  YY <- NULL
  lat <- NULL
  lon <- NULL
  alt <- NULL
  xy <- list(x=NULL,y=NULL)
  unit <- NULL
  country <- NULL
  descr <- NULL
  no.find <- TRUE
}
getdnmi<-list(val=x,station=station,yy=YY,
                lat=lat,lon=lon,alt=alt,
                x.0E65N=xy$x,y.0E65N=xy$y,
                location=location,
                wmo.no=NA,
                start=NA,yy0=NA,ele=NA,
                obs.name=descr, unit=unit,country=country,
                quality=NA,found=!no.find,
                ref=paste('The Norwegian Meteorological Institute',
                  'climatological archive'))
class(getdnmi) <- c("station","monthly.station.record")
getdnmi
}
