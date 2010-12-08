getnarp <- function(stnr=NULL,location=NULL,lon=NULL,lat=NULL,stations=NULL,
                    silent=TRUE,ele=101) {
	  	

data(narp.meta)

# Element 101, Mean monthly air-temperature.
# Element 111, Mean maximum monthly air-temperature.
# Element 112, Absolute maximum monthly air-temperature.
# Element 121, Mean minimum monthly air-temperature.
# Element 122, Aboslute minimum air-temperature.
# Element 401, Mean monthly sea level pressure.
# Element 601, Mean monthly precipitation sum.
# Element 602, Highest monthly 1-day precipitation.
# Element 701, Mean monthly days with snow cover > 50 %
# Element 801, Mean monthly cloud cover.

v.name<-switch(as.character(ele[1]),
                     '101'='monthly mean T(2m)','111'='mean maximum monthly air-temperature.',
                     '112'='absolute maximum monthly air-temperature',
                     '122'='aboslute minimum air-temperature',
                     '401'='monthly mean SLP',
                     '601'='monthly precipitation sum','602'='highest monthly 1-day precipitation',
                     '701'='mean monthly days with snow cover > 50 %','801'='monthly mean cloud cover')
unit<-switch(as.character(ele[1]),
                     '101'='degree Celsius','111'='degree Celsius','112'='degree Celsius',
                     '122'='degree Celsius','401'='hPa',
                     '601'='mm','701'='days','801'='%')

  if (is.character(stnr) & is.null(location)) {location <- stnr; stnr <- NULL} else
      if (is.character(stnr)) stnr <- NULL
  if (is.numeric(location) & is.null(stnr)) {stnr <- location; location <- NULL} else
      if (is.numeric(location)) location <- NULL

  #if (!silent) print("Retrieving the data from URL http://projects.met.no/~narp/")
  #if (!silent) print("Retrieving the data from ~/data/narp/")
  #if (!silent) print("Please be patient")
  if (!is.null(stnr)) {
     locmatch <- is.element(narp.meta$stnr,stnr)
     if (sum(locmatch)>0) {
       location <- narp.meta$names[locmatch]
       if (!silent) print(paste("Found",narp.meta$names[locmatch],"stnr=",stnr," lon=",
                          narp.meta$lons[locmatch]," lat=",narp.meta$lats[locmatch],
                          " country=",narp.meta$countries[locmatch]))
     } else {
       print(paste("Sorry, did NOT find station number '",stnr,"'",sep=""))
       print(substr(lower.case(stations$location),2,nchar(location)+1))
     }
     if (length(stnr)>1) {
       print(paste("getnarp: Found dublicates:",narp.meta$names[locmatch],"stnr=",stnr," lon=",
                   narp.meta$lons[locmatch]," lat=",narp.meta$lats[locmatch],
                   " country=",narp.meta$countries[locmatch]))
       i <- as.numeric(readline(paste("Which of these ( 1 -",length(stnr),")? ")))
       stnr <- stnr[i]
     } 
  }  else if (!is.null(location)) {
     locmatch <- is.element(lower.case(substr(narp.meta$names,1,4)),
                            lower.case(substr(location,1,4)))
     if (sum(locmatch)>0) {
       stnr <- narp.meta$stnr[locmatch]
       if (!silent) print(paste("Found",narp.meta$names[locmatch],"stnr=",stnr," lon=",
                          narp.meta$lons[locmatch]," lat=",narp.meta$lats[locmatch],
                          " country=",narp.meta$countries[locmatch]))
     } else {
       print(paste("Did not find '",location,"'",sep=""))
       print(substr(lower.case(stations$location),2,nchar(location)+1))
     }
     if (length(stnr)>1) {
       print(paste("Found dublicates:",narp.meta$names[locmatch],"stnr=",stnr," lon=",
                   narp.meta$lons[locmatch]," lat=",narp.meta$lats[locmatch],
                   " country=",narp.meta$countries[locmatch]))
       i <- as.numeric(readline(paste("Which of these ( 1 -",length(stnr),")? ")))
       stnr <- stnr[i]
     } 
  }

  if (!is.null(lon) & !is.null(lat) & is.null(stnr)) {
     if (!silent) print(paste("Find the nearest station to ",lon,
                              "E and ",lat,"N.",sep=""))
     dist <- distAB(lon,lat,narp.meta$lons,narp.meta$lats)
     distmatch <- dist == min(dist,na.rm=TRUE)
     distmatch[is.na(distmatch)] <- FALSE
     if (sum(distmatch,na.rm=TRUE)>0) {
       stnr <- narp.meta$stnr[distmatch]
       locmatch <- distmatch
       if (!silent) print(paste("Found",narp.meta$names[locmatch],"stnr=",stnr," lon=",
                          narp.meta$lons[locmatch]," lat=",narp.meta$lats[locmatch],
                          " country=",narp.meta$countries[locmatch]))
     }
  }

  #print(stnr)
  if (!is.null(stnr)) {
    #fname <- paste("http://projects.met.no/~narp/narp",ele,".txt",sep="")
    data(narp)
    #fname <- paste("~/data/narp/narp",ele,".txt",sep="")
    #print(fname)
    #narp <- read.fwf(fname,widths=c(5,3,4,rep(5,12)))
    iii <- is.element(narp[,1],stnr) & is.element(narp[,2],ele)
    narp[narp <= -9990] <- NA
    #print(table(narp$V1))   
    #print(sum(iii))                      
    obs <- station.obj(x=narp[iii,4:15]/10,yy=narp[iii,3],
                       obs.name=v.name,unit=unit,
                       station=as.numeric(stnr),ele=ele,wmo.no=narp.meta$WMO.number[locmatch],
                       location=narp.meta$names[locmatch],lon=narp.meta$lons[locmatch],
                       lat=narp.meta$lats[locmatch],alt=NA,
                       country=narp.meta$countries[locmatch],
                       ref=paste("URL: 'http://thule.oulu.fi/narp/' and ",
                                 "'http://projects.met.no/~narp/data_index.html'"))
  
} else obs <- list(name=narp.meta$names,lon=narp.meta$lons,lat=narp.meta$lats,
                   countries=narp.meta$countries,number=narp.meta$stnr)
  invisible(obs)
}
