rm(list=ls())
library(clim.pact)
source("ds_one.R")

elems <- c(101,601)
scens <- c("sresa1b","sresb1","sresa2")
test <- FALSE

options(device="png")

rcm.locs <- c("oslo","berg","trom","stoc","hels","koeb","reyk","tors")


for (ele in elems) {
  stations1 <- avail.locs(ele=ele)$name[is.element(avail.locs(ele=ele)$ident,"NORDKLIM")]
  countries1  <- avail.locs(ele=ele)$country[is.element(avail.locs(ele=ele)$ident,"NORDKLIM")]
  countries1 <- strip(countries1)
# number 31 is 'Ship' and causes an error...
  stations1 <- stations1[-c(31)]
  countries1 <- countries1[-c(31)]
  srt.1 <- c(23,5,6,54,18,22,28,26,27,31,34,35,36,33)
  srt.1 <- c(srt.1,(1:length(stations1))[-srt.1])
  stations2 <- avail.locs(ele=ele)$name[is.element(avail.locs(ele=ele)$ident,"NACD")]
  stations3 <- getnarp()$number; stations3 <- stations3[-2]

  for (scen in scens) {
    is <- 0
    print("NORDKLIM")
    for (station in stations1) {
      is <- is + 1
      if (countries1[is]=="N") predictand <- "nordklim+metno" else
                               predictand <- "nordklim"
      do.rcm <- (1:length(rcm.locs))[is.element(lower.case(substr(station,1,4)),rcm.locs)]
      #print(do.rcm)
      if ((ele!=101) | (is.null(do.rcm )) | (length(do.rcm)==0)) do.rcm <- 0
      ds.one(ele=ele,scen=scen,predictand=predictand,station=station,
                         do.rcm=FALSE,test=test,silent=TRUE) 
    
   }
    print("NACD")
    for (station in stations2) {
      do.rcm <- (1:length(rcm.locs))[is.element(lower.case(substr(station,1,4)),rcm.locs)]
      if ((ele!=101) | (is.null(do.rcm )) | (length(do.rcm)==0)) do.rcm <- 0
      ds.one(ele=ele,scen=scen,predictand="nacd",station=station,do.rcm=FALSE) 
   }
    print("NARP")
    for (station in stations3) {
      ds.one(ele=ele,scen=scen,predictand="narp",station=station) 
   }


 }
}

source("ds_tomas.R")
