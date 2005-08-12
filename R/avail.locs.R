# R.E. Benestad, met.no, Oslo, Norway 22.05.2002
# rasmus.benestad@met.no
#-------------------------------------------------------------------
# Selection of NACD and NORDKLIMstations.

avail.locs <- function(ele=101) {

#  source("strip.R")

  nacd.meta<-read.table('data/appendix.2')
  nordklim.meta<-read.fwf( 'data/nordklim_station_catalogue_v1_0.prn',
                 skip=1,as.is=TRUE,fill=TRUE,
                 width=c(2,30,12,11,11,4,3,2,4,3,2,
                         9,rep(6,23)),
                  col.names=c("i","location","height","country",
                              "number","Lat.deg","Lat.min","N.S",
                              "Lon.deg","Lon.min","E.W",
                              "ele101","ele101E","ele111","ele111E","ele112","ele112E",
                              "ele113","ele113E","ele121","ele121E","ele122","ele122E",
                              "ele123","ele123E","ele401","ele401E","ele601","ele601E",
                              "ele602","ele602E","ele701","ele701E","ele801","ele801E"))
  nacd <- length(nacd.meta$V5[is.element(nacd.meta$V14,ele)])
  narp <- getnarp()
  nnarp <- length(narp$lon)
  iele <- eval(parse(text=paste("!is.na(nordklim.meta$ele",ele,")",sep="")))
  nnordklim <- sum(iele)
  loc.list <- c(as.character(nacd.meta$V5[is.element(nacd.meta$V14,ele)]),
                     strip(as.character(nordklim.meta$location[iele])),narp$name)
  lat.list <- c(nacd.meta$V6[is.element(nacd.meta$V14,ele)] +
                1/60*nacd.meta$V7[is.element(nacd.meta$V14,ele)],
                nordklim.meta$Lat.deg[iele] + 1/60 * nordklim.meta$Lat.min[iele],narp$lat)
  lon.list <- c(nacd.meta$V9[is.element(nacd.meta$V14,ele)] +
                1/60*nacd.meta$V10[is.element(nacd.meta$V14,ele)],
                nordklim.meta$Lon.deg[iele] + 1/60 * nordklim.meta$Lon.min[iele],narp$lon)
  ew.list <- c(abbreviate(nacd.meta$V11[is.element(nacd.meta$V14,ele)]),
               abbreviate(nordklim.meta$E.W[iele]),rep("E",nnarp))
  lon.list[ew.list=="W"] <- lon.list[ew.list=="W"] * -1
  con.list <- c(as.character(nacd.meta$V3[is.element(nacd.meta$V14,ele)]),
                as.character(nordklim.meta$country[iele]),narp$countries)
  avail.locs<-list(name=loc.list,
                   lons=lon.list,
                   lats=lat.list,
                   country=factor(strip(abbreviate(con.list))),
                   nacd=nacd,nnordklim=nnordklim,
                   ident=c(rep("NACD",nacd),rep("NORDKLIM",nnordklim),rep("NARP",nnarp))) 
  avail.locs
}
