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


stationmap <- function(ele=101,NORDKLIM=TRUE,NACD=TRUE,silent=TRUE) {

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

plot(c(-80,40),c(50,82),type="n",
     main=ele.c,xlab="Longitude",ylab="Latitude")
addland()
if (!silent) print(locs)

for (loc in locs) {
#  print("NACD:")
  if ((NORDKLIM) & (NACD)) { 
    obs.nacd<-getnacd(loc,silent=TRUE)
    if (obs.nacd$found) {
      points(obs.nacd$lon,obs.nacd$lat,col="blue",pch=20,cex=1.25)
    }
#     print("NordKlim")
    obs.nork<-getnordklim(loc,silent=TRUE)
    if (obs.nork$found) {
        points(obs.nork$lon,obs.nork$lat,col="red",pch=20,cex=0.8)
    }
  } else if (NORDKLIM) {
    obs.nork<-getnordklim(loc,silent=TRUE)
      if (obs.nork$found) {
        points(obs.nork$lon,obs.nork$lat,col="red",pch=20,cex=0.8)
      }
  } else if (NACD) {
    obs.nacd<-getnacd(loc,silent=TRUE)
    if (obs.nacd$found) {
      points(obs.nacd$lon,obs.nacd$lat,col="blue",pch=20,cex=1.25)
    }
  }
}
grid()
legend(-75,55,c('NACD','NordKlim'),pch=c(20,20),
       col=c('blue','red'),bg="grey95")

}
