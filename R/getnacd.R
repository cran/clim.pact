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

getnacd <- function(location="prompt",ele.c='101',ascii=FALSE,silent=FALSE,
                    direc="data") {

if (location=="prompt") {
  locs <- avail.locs(as.integer(ele.c))
  print(length(locs))
  locs.name <- locs$name[locs$ident=="NACD"]
  print(locs.name)
  i.loc <- readline(prompt="Please enter the number of desired location: ")
  i.loc <- as.integer(i.loc)
  location <- locs.name[i.loc]
} 


if (is.character(ele.c)) {
  ele.c<-lower.case(ele.c)
  ele.c<-switch(ele.c,
                't2m'='101','rr'='601','slp'='401','cloud'='801',
                't2'='101','precip'='601','101'='101','401'='401',
                '601'='601','801'='801')
} else {
  ele.c<-as.character(ele.c)
}

#print(location)
#print(ele.c)
fr.name<-paste(direc,'/getnacd_',ele.c,'.Rdata',sep="")
ascii<- ascii | !file.exists(fr.name)

if (ascii) {
# Read the original ASCII files - slow
  obs<-read.fwf(
     paste('data/nacd_v1.',ele.c,sep=""),
     width=c(5,3,4,rep(5,12),3))

# Read the information about the stations: Metadata

  data(meta.nacd)

# Save as R-data-file
  save(obs,meta.nacd,file=fr.name)
}

# Load R-data-file - FAST!

load(fr.name)

station<-obs$V1
ele<-obs$V2

yy<-obs$V3
country<-obs$V16

scale<-10
if (ele[1]==801) scale<-1

val<-as.matrix(obs[,4:15])/scale
val[val <= -99.9] <- NA

#print(obs[1,])
#print(meta.nacd[1,])
#print(c(as.character(meta.nacd$V1[1]),as.character(meta.nacd$V2[1]),
#        as.character(meta.nacd$V3[1]),as.character(meta.nacd$V4[1]),
#        as.character(meta.nacd$V5[1]),as.character(meta.nacd$V6[1]),
#        as.character(meta.nacd$V7[1]),as.character(meta.nacd$V8[1]),
#        as.character(meta.nacd$V9[1]),as.character(meta.nacd$V10[1]),
#        as.character(meta.nacd$V11[1]),as.character(meta.nacd$V12[1]),
#        as.character(meta.nacd$V13[1]),as.character(meta.nacd$V14[1]),
#        as.character(meta.nacd$V15[1]),as.character(meta.nacd$V16[1])))

if (is.character(location)) {
  nc<-nchar(location)
  location<-paste(upper.case(location),
                paste(rep(" ",21-nc),sep="",collapse=""),sep="")

  no.find<-FALSE
  if ((sum(is.element(meta.nacd$location,location) &
            is.element(as.numeric(as.character(meta.nacd$element)),ele))==0) &
      !(silent)) no.find<-TRUE
  meta<-meta.nacd[is.element(meta.nacd$location,location) &
                  is.element(as.numeric(as.character(meta.nacd$element)),ele),]

} else if (is.numeric(location)){
  if ((sum(is.element(meta.nacd$station.number,location) &
            is.element(as.numeric(as.character(meta.nacd$element)),ele))==0) &
      !(silent)) no.find<-TRUE
  meta<-meta.nacd[is.element(meta.nacd$station.number,location) &
                  is.element(as.numeric(as.character(meta.nacd$element)),ele),]
}

if (no.find) {
  print("getnacd: ERROR - cannot find the right record!")
  print(sum(is.element(meta.nacd$location,location) &
              is.element(meta.nacd$element,ele)))

  print("location")
  print(location)
  print("levels(meta.nacd$location)")
  print(levels(meta.nacd$location))

  print("ele")
  print(table(ele))
  print("levels(meta.nacd$element)")
  print(table(as.numeric(as.character(meta.nacd$element))))

  print(paste("sum(is.element(meta.nacd$location,location))=",
              sum(is.element(meta.nacd$location,location))))
  print(paste("sum(is.element(meta.nacd$element,ele))=",
              sum(is.element(meta.nacd$element,ele)))) 
  print("meta:")
  print(meta)
  print("station:")
  print(summary(station))

  print("country:")
  print(levels(country))
  print(meta$country)
}
 
iloc<-is.element(station,meta$station.number) &
                (country == meta$country) 

if (sum(iloc)==0) {
 print("summary(iloc)")
 print(summary(iloc))
 print("sum(iloc)")
 print(sum(iloc))
 print("sum(is.element(station,meta$station.number))")
 print(sum(is.element(station,meta$station.number)))
 print("sum(country == meta$country)")
 print(sum(country == meta$country))
}

obs.name<-switch(as.character(ele[1]),
                     '101'='monthly mean T(2m)',
                     '401'='monthly mean SLP',
                     '601'='monthly precipitation sum',
                     '801'='monthly mean cloud cover')
unit<-switch(as.character(ele[1]),
                     '101'='degree Celsius',
                     '401'='hPa',
                     '601'='mm',
                     '801'='%')
#print(as.character(meta$V16))
quality<-switch(as.character(meta$quality),
                ' H'='Homogenous, rigorously tested & adjusted',
                'H'='Homogenous, rigorously tested & adjusted',
                ' T'='Tested, maybe adjusted but not perfectly H.',
                'T'='Tested, maybe adjusted but not perfectly H.',
                ' N'='Not tested for inhomogenouity',
                'N'='Not tested for inhomogenouity',
                ' E'='Environm. changes prevents clim.change studies',
                'E'='Environm. changes prevents clim.change studies',
                ' I'='Inhomogenous series which presently are unadjustable',
                'I'='Inhomogenous series which presently are unadjustable')

lat<-meta$degN + meta$minN/60
lon<-meta$degE + meta$minE/60
lat[meta$N.S==" S"]<-lat[meta$N.S==" S"]*-1
lon[meta$E.W==" W"]<-lon[meta$E.W==" W"]*-1
#print(levels(meta$V8))
#print(levels(meta$V11))

xy<-COn0E65N(lon,lat)

getnacd<-list(val=val[iloc,],station=meta$station.number,yy=yy[iloc],
              lat=lat,lon=lon,alt=meta$alt,
              x.0E65N=xy$x,y.0E65N=xy$y,
              location=location, wmo.no=meta$wmo.number,
              start=meta$start,yy0=meta$year.1,ele=ele[1],
              obs.name=obs.name, unit=unit,country=meta$country,
              quality=quality,found=!no.find,
              ref='Frich et al. (1996), DMI scientific report 96-1')
class(getnacd) <- c("station","monthly.station.record")
getnacd
}
