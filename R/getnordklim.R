# This R routine reads the NACD data. The code
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

getnordklim <- function(location="prompt",ele.c='101',
                        ascii=FALSE,silent=FALSE,direc="data") {


#library(stringfun)
#source("COn0E65N.R")
#source("avail.locs.R")

#if (location=="prompt") {
#  locs <- avail.locs(as.integer(ele.c))
#  print(length(locs))
#  locs.name <- locs$name[locs$ident=="NORDKLIM"]
#  print(locs.name)
#  i.loc <- readline(prompt="Please enter the number of desired location: ")
#  i.loc <- as.integer(i.loc)
#  location <- locs.name[i.loc]
#}
ele <- c(101,111,112,113,121,122,123,401,601,602,701,801,911)

#print(ele.c)  
if (is.character(ele.c)) {
  ele.c<-lower.case(ele.c)
  ele.c<-switch(ele.c,
                't2m'='101','rr'='601','slp'='401','cloud'='801',
                't2'='101','precip'='601','101'='101','401'='401','111'='111',
                '112'='112','113'='113','121'='121','122'='122','123'='123',
                '601'='601','602'='602','801'='801','911'='911')
} else {
  ele.c<-as.character(ele.c)
}

#print(location)
#print(ele.c)
fr.name<-paste(direc,'/nordklim_',ele.c,'.Rdata',sep="")
ascii<- ascii | !file.exists(fr.name)

if (ascii) {
# Read the original ASCII files - slow
  obs<-read.table(paste('data/p',ele.c,'.prn',sep=""))

# Read the information about the stations: Metadata
# 

  data(nordklim.meta)

# Save as R-data-file
  save(obs,meta,file=fr.name)
}

# Load R-data-file - FAST!

load(fr.name)

station<-as.integer(obs$V1)
ele<-obs$V2
yy<-obs$V3
country<-strip(obs$V16)

obs.name<-switch(as.character(ele[1]),
              '101'='mean T(2m)',
              '111'='mean maximum T(2m)',
              '112'='highest maximum T(2m)',
              '113'='day of Th date Thd',
              '121'='mean minimum T(2m)',
              '122'='lowest minimum T(2m)',
              '123'='day of Tl date Tld',
              '401'='mean SLP',
              '601'='monthly accum. precip.',
              '602'='maximum precip.',
              '701'='Number of days with snow cover (> 50% covered) days dsc',
              '801'='Mean cloud cover % N',
              '911'='mean snow depth')

unit <-switch(as.character(ele[1]),
              '101'='deg C',
              '111'='deg C',
              '112'='deg C',
              '113'='date',
              '121'='deg C',
              '122'='deg C',
              '123'='date',
              '401'='hPa',
              '601'='mm/month',
              '602'='mm/day',
              '701'='days',
              '801'='%',
              '911'='cm')
#print(paste("obs name=",obs.name,"; unit=",unit))

scale <- switch(as.character(ele[1]),
               '101'=0.1,
               '111'=0.1,
               '112'=0.1,
               '113'=1,
               '121'=0.1,
               '122'=0.1,
               '123'=1,
               '401'=0.1,
               '601'=0.1,
               '602'=0.1,
               '701'=1,
               '801'=1,
               '911'=1)
#print("scale")
#print(scale)
#print(ele[1])
#print(dim(obs))
#print(summary(obs))

val<-as.matrix(obs[,4:15])*scale
val[val <= -99.9] <- NA

#print(paste("GETNORDKLIM: ",location))
location<-upper.case(substr(location,1,4))

#print("search appendix:")
#print(as.character(meta$V5))
locations <- upper.case(substr(as.character(meta$location),1,4))
in.app<-is.element(locations,location) 
#print("in.app=")
#print(sum(in.app))
location <- strip(as.character(meta$location[in.app]))

no.find<-FALSE
if ((sum(in.app)==0) & !(silent)) {
  print("getnordklim: ERROR - cannot find the right record!")
  print(sum(is.element(meta$location,location)))
  
  print("ele.c")
  print(ele.c)
  print("location")
  print(location)
  print("table(ele)")
  print(table(ele))
  print(as.character(meta$location))
  print("sum(is.element(meta$location,location))")
  print(sum(is.element(meta$location,location)))
  print("meta:")
  print(meta)
  print("station:")
  print(summary(station))

  print("country:")
  print(levels(country))
  print(meta$V3)
  no.find<-TRUE

  print("Available locations:")
  print(meta$location)

  
#if (no.find) {
# print("summary(iloc)")
# print(summary(iloc))
# print("sum(iloc)")
# print(sum(iloc))
# print("sum(is.element(station,meta$V2))")
# print(sum(is.element(station,meta$V2)))
# print("sum(country == meta$V3)")
# print(sum(country == meta$V3))
#}
}

#print(length(!no.find))
if (!no.find) {
  meta<-meta[in.app,]
#print(meta)
#print(table(station))
#print(meta$number)

  iloc<-is.element(station,meta$number)

#print(sum(iloc))
#print(as.character(meta$V16))

#  print(" Test 1:")
#  print(meta[1,])
#  print(c(meta$Lat.deg,meta$Lat.min,meta$Lon.deg,meta$Lon.min))

  lat<-meta$Lat.deg + meta$Lat.min/60
  lon<-meta$Lon.deg + meta$Lon.min/60
#  print(paste("GETNORDKLIM: call strip for N.S & E.W","'",
#              meta$N.S,"', '",meta$E.W,"'",sep=""))
  meta$N.S<-strip(meta$N.S)
  meta$E.W<-strip(meta$E.W)
  lat[(meta$N.S=="S") | (meta$N.S==" S")]<-lat[(meta$N.S=="S") | (meta$N.S==" S")]*-1
  lon[(meta$E.W=="W") | (meta$E.W==" W")]<-lon[(meta$E.W=="W") | (meta$E.W==" W")]*-1
  
#print(levels(meta$V8))
#print(levels(meta$V11))

  xy<-COn0E65N(lon,lat)
} else {
  lat<-NA
  lon<-NA
  xy<-list(x=NA,y=NA)
  location<-"Not found"
}
#print(" Test 2: dim(val)")
#print(length(station))
#print(length(yy))
#print(dim(val))
#print(dim(val[iloc,]))
#print(length(yy[iloc]))
#print(" Test 3:")
#print(c(lon,lat,meta$altitude,ele[1]))

getnordklim<-list(val=val[iloc,],station=meta$number,yy=yy[iloc],
              lat=lat,lon=lon,alt=meta$height,
              x.0E65N=xy$x,y.0E65N=xy$y,
              location=location, wmo.no=NA,
              start=NA,yy0=NA,ele=ele[1],
              obs.name=obs.name, unit=unit,country=meta$country,
              found=!no.find,
              ref=paste('Rissanen et al., (2000), DNMI KLIMA 10/00,',
                        'Norwegian Meteololog. Inst., met.no'))
class(getnordklim) <- c("station","monthly.station.record")
getnordklim
}
