rm(list=ls())
library(clim.pact)
source("clim.pact/R/ds.R")

ele <- 101
cmon <- "Jan"
scen <- "_B2"

eof.list <- avail.eofs()
eof.list <- eof.list[grep(cmon,eof.list)]
eof.list <- eof.list[grep(scen,eof.list)]
eof.list <- eof.list[grep("ncep",eof.list)]
eof.list <- eof.list[grep("_Mon.Rdata",eof.list)]
print(eof.list)
station.list <- avail.locs(ele)
station.list$name <- station.list$name[station.list$name=="OSLO-BLINDERN"]
station.list$name <- station.list$name[!is.na(station.list$name)]
print(station.list$name)
for (i.loc in 1:length(station.list$name)) {
  location <- station.list$nam[i.loc]
  dataset <- station.list$ident[i.loc]
  if (dataset=="NACD") obs<-getnacd(location,ele) else
  if (dataset=="NORDKLIM") obs<-getnordklim(location,ele)

  if (!is.null(range(obs$yy))) {
    for (i.eof in 1:length(eof.list)) {
      load(paste("data/",eof.list[i.eof],sep=""))
      ds(preds=eof,dat=obs,plot=F)
      rm(eof)
    }
  }
}
