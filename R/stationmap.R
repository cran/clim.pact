# Empirical downscaling using "mixed common EOFs" from ceof.R
# Predictand is a time series from NACD or climate station.
# Monthly mean values.
#
# Reference: R.E. Benestad et al. (2002),
#            Empirically downscaled temperature scenarios for Svalbard,
#            submitted to Atm. Sci. Lett.
#
#            R.E. Benestad (2001),
#            A comparison between two empirical downscaling strategies,
#            Int. J. Climatology, 1645-1668, vol. 21, DOI 10.1002/joc.703
#
# R.E. Benestad, met.no, Oslo, Norway 16.04.2002
# rasmus.benestad@met.no
#------------------------------------------------------------------------


stationmap <- function(ele=101,NORDKLIM=TRUE,NACD=TRUE,silent=TRUE,names=FALSE,
                       name.len=4,x.offs=0.1,y.offs=-0.5,str.cex=0.7,
                       countries=NULL,x.rng=NULL,y.rng=NULL) {

# Load libraries, and compile function:
  
#source("getnordklim.R")
#source("getnacd.R")
#source("addland.R")
#source("strip.R")

ele.c<-switch(as.character(ele),
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
              '901'='mean snow depth')

#nacd.meta<-read.table('data/appendix.2')
#nordklim.meta<-read.fwf( 'data/nordklim_station_catalogue_v1_0.prn',
#                        skip=1,width=c(2,30,12,11,11,4,3,2,4,3,2),
#                        col.names=c("i","location","height","country",
#                                    "number","Lat.deg","Lat.min","N.S",
#                                    "Lon.deg","Lon.min","E.W"))
data(nacd.meta)
data(nordklim.meta)

#-------------------------------------------------------------------
# Selection of NACD and NORDKLIMstations.

if ((NORDKLIM) & (NACD)) { 
 locs<-c(as.character(nacd.meta$V5[is.element(nacd.meta$V14,ele)]),
         strip(as.character(meta$location)))
} else if (NORDKLIM) locs<-as.character(strip(as.character(meta$location))) else
       if (NACD) locs<-as.character(nacd.meta$V5[is.element(nacd.meta$V14,ele)])

keep <- rep(TRUE,length(locs))
if (!is.null(countries)) {
  countries <- switch(lower.case(countries),
                    "sweden"="S","sverige"="S","finland"="FIN","denmark"="DK",
                    "danmark"="DK","norway"="N","norge"="N","noreg"="N",
                    "island"="IS","iceland"="IS","faeroe islands"="FR",
                    "belgium"="B","greenland"="G","great britain"="GB",
                    "united kingdom"="GB","uk"="GB","u.k."="GB",
                    "ireland"="IRL","netherlands"="NL","holland"="NL")
  keep <- keep & c(is.element(strip(nacd.meta$V3),countries),
                   is.element(strip(meta$country),countries))
}

locs <- locs[keep]

plot(c(-80,40),c(50,82),type="n",
     main=ele.c,xlab="Longitude",ylab="Latitude")
addland()
if (!silent) print(locs)

for (loc in locs) {
  if ((NORDKLIM) & (NACD)) { 
    obs.nacd<-getnacd(loc,silent=TRUE)
    if (obs.nacd$found) {
      points(obs.nacd$lon,obs.nacd$lat,col="blue",pch=20,cex=1.25)
      if (names) text(obs.nacd$lon+x.offs,obs.nacd$lat+y.offs,col="blue",
                      substr(loc,1,name.len),cex=str.cex,pos=4)
    }
    obs.nork<-getnordklim(loc,silent=TRUE)
    if (obs.nork$found) {
      points(obs.nork$lon,obs.nork$lat,col="red",pch=20,cex=0.8)
      if (names) text(obs.nacd$lon+x.offs,obs.nacd$lat+y.offs,col="red",
                      substr(loc,1,name.len),cex=str.cex,pos=4)
    }
  } else if (NORDKLIM) {
    obs.nork<-getnordklim(loc,silent=TRUE)
      if (obs.nork$found) {
      points(obs.nork$lon,obs.nork$lat,col="red",pch=20,cex=0.8)
      if (names) text(obs.nacd$lon+x.offs,obs.nacd$lat+y.offs,col="red",
                      substr(loc,1,name.len),cex=str.cex,pos=4)
      }
  } else if (NACD) {
    obs.nacd<-getnacd(loc,silent=TRUE)
    if (obs.nacd$found) {
      points(obs.nacd$lon,obs.nacd$lat,col="blue",pch=20,cex=1.25)
      if (names) text(obs.nacd$lon+x.offs,obs.nacd$lat+y.offs,col="blue",
                      substr(loc,1,name.len),cex=str.cex,pos=4)
    }
  }
}
grid()
legend(-75,55,c('NACD','NordKlim'),pch=c(20,20),
       col=c('blue','red'),bg="grey95")

}
