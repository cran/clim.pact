getnarp <- function(stnr=NULL,location=NULL,lon=NULL,lat=NULL,stations=NULL,silent=TRUE,ele=101) {
	  	

narp.names <- c("Torshavn","Strond Kr.st.","Upernavik","Ilulissat Airport","Nuuk","Narsarsuaq",
                "Danmarkshavn","Ittoqqortoormiit","Tasiilaq","Stykkisholmur","Reykjavik","Vestmannaeyar",
                "Haell","Akureyri","Raufarhoefn","Teigarhorn","Ship M","Tromsoe","Vardoe","Bjoernoeya",
                "Hopen","Svalbard Airport","Svalbard Airport reconstructed","Jan Mayen")
narp.lons <- c( -6.7667, -6.5833,-56.1667,-51.0667,-51.7500,-48.1667,-18.7667,-22.0000,-37.6333,
               -22.7333,-21.9000,-20.2833,-20.4167,-18.0833,-15.9500,-15.2000,  2.0000, 18.9333,
                31.0833, 19.0167, 25.0667, 15.4667, 15.4667, -8.6667)
narp.lats <- c( 62.0167, 62.2667, 72.7833, 69.2330, 64.1667, 61.2000, 76.7667, 70.4833, 65.6000,
                65.0833, 64.1333, 63.4000, 64.0667, 65.6833, 66.4500, 64.3000, 66.0000, 69.6500,
                70.3667, 74.5167, 76.5000, 78.2500, 78.2500, 70.9333)
narp.countries <- c(rep("Faeroes",2),rep("Greenland",7),rep("Iceland",7),rep("Norway",8))
who.stnr <- c(6011,NA,4210,4221,4250,4270,4320,4339,4360,4013,4030,4048,4050,4063,4077,4092,NA,
              1026,1098,1028,1062,1008,NA,1001)
narp.stnr <- c(06011,33054,04210,04221,04250,04270,04320,04339,04360,04013,04030,04048,04050,
               04063,04077,04092,76900,90450,98550,99710,99720,99840,99841,99950)

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

  if (!silent) print("Retrieving the data from URL http://projects.met.no/~narp/")
  if (!silent) print("Please be patient")
  if (!is.null(stnr)) {
     locmatch <- is.element(narp.stnr,stnr)
     if (sum(locmatch)>0) {
       location <- narp.names[locmatch]
       if (!silent) print(paste("Found",narp.names[locmatch],"stnr=",stnr," lon=",
                          narp.lons[locmatch]," lat=",narp.lats[locmatch],
                          " country=",narp.countries[locmatch]))
     } else {
       print(paste("Sorry, did NOT find station number '",stnr,"'",sep=""))
       print(substr(lower.case(stations$location),2,nchar(location)+1))
     }
     if (length(stnr)>1) {
       print(paste("Found dublicates:",narp.names[locmatch],"stnr=",stnr," lon=",
                   narp.lons[locmatch]," lat=",narp.lats[locmatch],
                   " country=",narp.countries[locmatch]))
       i <- as.numeric(readline(paste("Which of these ( 1 -",length(stnr),")? ")))
       stnr <- stnr[i]
     } 
  }  else if (!is.null(location)) {
     locmatch <- is.element(lower.case(substr(narp.names,1,4)),lower.case(substr(location,1,4)))
     if (sum(locmatch)>0) {
       stnr <- narp.stnr[locmatch]
       if (!silent) print(paste("Found",narp.names[locmatch],"stnr=",stnr," lon=",
                          narp.lons[locmatch]," lat=",narp.lats[locmatch],
                          " country=",narp.countries[locmatch]))
     } else {
       print(paste("Did not find '",location,"'",sep=""))
       print(substr(lower.case(stations$location),2,nchar(location)+1))
     }
     if (length(stnr)>1) {
       print(paste("Found dublicates:",narp.names[locmatch],"stnr=",stnr," lon=",
                   narp.lons[locmatch]," lat=",narp.lats[locmatch],
                   " country=",narp.countries[locmatch]))
       i <- as.numeric(readline(paste("Which of these ( 1 -",length(stnr),")? ")))
       stnr <- stnr[i]
     } 
  }

  if (!is.null(lon) & !is.null(lat) & is.null(stnr)) {
     if (!silent) print(paste("Find the nearest station to ",lon,
                              "E and ",lat,"N.",sep=""))
     dist <- distAB(lon,lat,narp.lons,narp.lats)
     distmatch <- dist == min(dist,na.rm=TRUE)
     distmatch[is.na(distmatch)] <- FALSE
     if (sum(distmatch,na.rm=TRUE)>0) {
       stnr <- narp.stnr[distmatch]
       locmatch <- distmatch
       if (!silent) print(paste("Found",narp.names[locmatch],"stnr=",stnr," lon=",
                          narp.lons[locmatch]," lat=",narp.lats[locmatch],
                          " country=",narp.countries[locmatch]))
     }
  }

  #print(stnr)
  if (!is.null(stnr)) {
    fname <- paste("http://projects.met.no/~narp/narp",ele,".txt",sep="")
    print(fname)
    narp <- read.fwf(fname,widths=c(5,3,4,rep(5,12)))
    iii <- is.element(narp$V1,stnr)
    narp[narp <= -9990] <- 0
    #print(table(narp$V1))   
    #print(sum(iii))                      
    obs <- station.obj(x=as.matrix(narp[iii,4:15])/10,yy=narp$V3[iii],
                       obs.name=v.name,unit=unit,
                       station=as.numeric(stnr),ele=ele,,wmo.no=who.stnr[locmatch],
                       location=narp.names[locmatch],lon=narp.lons[locmatch],
                       narp.lats[locmatch],alt=NA,country=narp.countries[locmatch],
                       ref=paste("URL: 'http://thule.oulu.fi/narp/' and ",
                                  "'http://projects.met.no/~narp/data_index.html'"))
  
} else obs <- list(name=narp.names,lon=narp.lons,lat=narp.lats,countries=narp.countries,
                   number=narp.stnr)
  invisible(obs)
}
