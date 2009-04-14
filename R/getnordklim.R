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

getnordklim <- function(location=NULL,ele.c='101',
                        ascii=FALSE,silent=FALSE,direc="data") {


#library(stringfun)
#source("COn0E65N.R")
#source("avail.locs.R")

if (is.null(location)) {
  locs <- avail.locs(as.integer(ele.c))
  print(length(locs))
  locs.name <- locs$name[locs$ident=="NORDKLIM"]
  locs.name[is.element(upper.case(locs.name),"JAN")] <- "Jan_Mayen"
  return(locs.name)
}
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

  data(nordklim.meta,envir = environment())

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

val<-as.matrix(obs[,4:15])*scale
val[val <= -99.9] <- NA
#location<-upper.case(substr(location,1,4))
#locations <- upper.case(substr(as.character(meta$location),1,4))
#in.app<-is.element(locations,location)
test.char <- 4
Location <- location
in.app <- rep(F,length(meta$location))
while ( (sum(in.app)!=1) & (test.char <= nchar(Location)) ) {
  location<-upper.case(substr(Location,1,test.char))
  locations <- upper.case(substr(as.character(meta$location),1,test.char))
  in.app<-is.element(locations,location)
  test.char <- test.char +1
}
#print(c(sum(in.app)!=1,test.char <= nchar(Location)))
#stop(paste(sum(in.app),"-",test.char))

no.find<-FALSE
if ((sum(in.app)==0) & !(silent)) {
  print("getnordklim: ERROR - cannot find the right record!")
  print(sum(is.element(meta$location,location)))
  print("ele.c");  print(ele.c)
  print("location"); print(location)
  print("table(ele)"); print(table(ele))
  print(as.character(meta$location))
  print("sum(is.element(meta$location,location))")
  print(sum(is.element(meta$location,location)))
  print("meta:"); print(meta)
  print("station:"); print(summary(station))
  print("country:"); print(levels(country))
  print(meta$V3)
  no.find<-TRUE
  print("Available locations:")
  print(meta$location)
} else if (sum(in.app)==1) {
  location <- strip(as.character(meta$location[in.app]))
  meta<-meta[in.app,]
  iloc<-is.element(station,meta$number)
  lat<-meta$Lat.deg + meta$Lat.min/60
  lon<-meta$Lon.deg + meta$Lon.min/60
  meta$N.S<-strip(meta$N.S)
  meta$E.W<-strip(meta$E.W)
  lat[(meta$N.S=="S") | (meta$N.S==" S")]<-lat[(meta$N.S=="S") | (meta$N.S==" S")]*-1
  lon[(meta$E.W=="W") | (meta$E.W==" W")]<-lon[(meta$E.W=="W") | (meta$E.W==" W")]*-1

  xy<-COn0E65N(lon,lat)

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
} else {
  getnordklim<-list(location=location,found=FALSE)
}
getnordklim
}
