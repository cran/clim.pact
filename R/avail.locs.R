# R.E. Benestad, met.no, Oslo, Norway 22.05.2002
# rasmus.benestad@met.no
#-------------------------------------------------------------------
# Selection of NACD and NORDKLIMstations.

avail.locs <- function(ele=101) {

#  source("strip.R")

  nacd.meta<-read.table('data/appendix.2')
  nordklim.meta<-read.fwf( 'data/nordklim_station_catalogue_v1_0.prn',skip=1,
                 width=c(2,30,12,11,11,4,3,2,4,3,2),
                  col.names=c("i","location","height","country",
                              "number","Lat.deg","Lat.min","N.S",
                              "Lon.deg","Lon.min","E.W"))
  nacd <- length(nacd.meta$V5[is.element(nacd.meta$V14,ele)])
  nnordklim <- length(nordklim.meta$location)
  loc.list <- c(as.character(nacd.meta$V5[is.element(nacd.meta$V14,ele)]),
                     strip(as.character(nordklim.meta$location)))
  lat.list <- c(nacd.meta$V6[is.element(nacd.meta$V14,ele)] +
                1/60*nacd.meta$V7[is.element(nacd.meta$V14,ele)],
                nordklim.meta$Lat.deg + 1/60 * nordklim.meta$Lat.min)
  lon.list <- c(nacd.meta$V9[is.element(nacd.meta$V14,ele)] +
                1/60*nacd.meta$V10[is.element(nacd.meta$V14,ele)],
                nordklim.meta$Lon.deg + 1/60 * nordklim.meta$Lon.min)
  ew.list <- c(abbreviate(nacd.meta$V11[is.element(nacd.meta$V14,ele)]),
               abbreviate(nordklim.meta$E.W))
  lon.list[ew.list=="W"] <- lon.list[ew.list=="W"] * -1
  con.list <- c(as.character(nacd.meta$V3[is.element(nacd.meta$V14,ele)]),
                as.character(nordklim.meta$country))
  avail.locs<-list(name=loc.list,
                   lons=lon.list,
                   lats=lat.list,
                   country=factor(strip(abbreviate(con.list))),
                   nacd=nacd,nnordklim=nnordklim,
                   ident=c(rep("NACD",nacd),rep("NORDKLIM",nnordklim))) 
  avail.locs
}
