getClimateExplorer <- function(location=NULL,URL="http://climexp.knmi.nl/data/",ele.c="t") {
  if (is.null(location)) {
    element <- switch(ele.c,"t"="temperature",
                      "p"="precipitation","l"="sealev",)
    list <-  readLines(paste("http://climexp.knmi.nl/allstations.cgi?someone@somewhere+",
                                element,"+12",sep=""))
    iwmo <- grep("WMO station code:",list)
    icrd <- grep("coordinates:",list)
    iyrs <- grep("years with data in",list)
    iter <- grep("Terrain:",list)

    wmo <- as.numeric(substr(list[iwmo],18,24))
    coord <- substr(list[icrd],18,nchar(list[icrd]))
    name <- list[icrd-1]
    nyrs <- as.numeric(substr(list[iyrs],6,10))
    situation <- list[iwmo+1]
    N <- length(wmo)
    country <- rep("",N); lon <- rep(NA,N); lat <- lon; alt <- lon; pop <- lon
    for (i in 1:N) {
      paren1 <- instring("(",name[i]); paren2 <- instring(")",name[i]);
      country[i] <- substr(name[i],paren1+1,paren2-1)
      non.chars <- c( instring("(",name[i]),instring(",",name[i]),
              instring("   ",name[i]) )
      non.chars <- non.chars[!non.chars==0]
      eon <- min(non.chars)
      if (eon[1]>0) location <- substr(name[i],1,eon[1]-1) else
                    location <- strip(substr(name[i],1,nchar(name[i])))
      while (substr(location,nchar(location),nchar(location))==" ")
        location <- substr(location,1,nchar(location)-1)    
      if (substr(location,1,1)=="-") location <- substr(location,2,nchar(location))
      spaces <- instring(" ",location)
      if (length(spaces)>0) {
        for (ij in 1:length(spaces))
          if (spaces[ij]>1) location <- paste(substr(location,1,spaces[ij]-1),
                                              "_",substr(location,spaces[ij]+1,
                                              nchar(location)),sep="") else
          if (spaces[ij]==1) location <- substr(location,2,nchar(location))
      }

      apos <- instring("'",location)
      if (length(apos)>0) {
         for (ij in 1:length(apos))
           location <- paste(substr(location,1,apos[ij]-1),
                          substr(location,apos[ij]+1,nchar(location)),sep="")
      }

      slash.hunt <- instring("/",location)
      if (length(slash.hunt)>0) {
        for (iii in 1:length(slash.hunt)) {
          if (slash.hunt[iii]>0) location <-
            paste(substr(location,1,slash.hunt[1]-1),
            substr(location,slash.hunt[1]+1,nchar(location)),sep="-")
        }
      }
      name[i] <- location
      i1n <- instring("N,",coord[i])[1]
      i1s <- instring("S,",coord[i])[1]
      i1 <- max(c(i1n,i1s),na.rm=TRUE)
      i2e <- instring("E,",coord[i])[1]
      i2w <- instring("W,",coord[i])[1]
      i2 <- max(c(i2e,i2w),na.rm=TRUE)
      i3 <-  instring("m ",coord[i])[1]
      lat[i] <- as.numeric(substr(coord[i],15,i1-1))
      lon[i] <- as.numeric(substr(coord[i],i1+3,i2-1))
      alt[i] <- as.numeric(substr(coord[i],i2+3,i3-1))
      paren1 <- instring("(",situation[i]); paren2 <- instring(")",situation[i]);
      if (paren1[1]==0)paren1 <- nchar(situation[i])
      if (paren2[1]==0)paren2 <- nchar(situation[i])
      br <- instring("<br>",situation[i])
      ipop <- instring("pop",situation[i])
      if (length(ipop)>0) pop[i] <- as.numeric(substr(situation[i],ipop+4,paren2))
      if (length(instring("Associated with ",situation[i]))>0)
        situation[i] <- substr(situation[i],18,paren1) else
        situation[i] <- substr(situation[i],18,br)
    }
    results <- data.frame(wmo=wmo,name=name,situation=situation,lon=lon,lat=lat,
                          alt=alt,n.years=nyrs,pop=pop)
    return(results)
  }

  fullURL <- paste(URL,ele.c,location,".dat",sep="")
  print(fullURL)
  header <- readLines(fullURL,n=5)
  print(header)

  # Meta-data:
  # Tidy up the name string:

      non.chars <- c( instring("(",header[2]),instring(",",header[2]),
              instring("   ",header[2]) )
      non.chars <- non.chars[!non.chars==0]
      eon <- min(non.chars)
      if (eon[1]>0) location <- substr(header[2],1,eon[1]-1) else
                    location <- strip(substr(header[2],1,nchar(header[2])))
      while (substr(location,nchar(location),nchar(location))==" ")
        location <- substr(location,1,nchar(location)-1)    
      if (substr(location,1,1)=="-") location <- substr(location,2,nchar(location))
      spaces <- instring(" ",location)
      if (length(spaces)>0) {
        for (ij in 1:length(spaces))
          if (spaces[ij]>1) location <- paste(substr(location,1,spaces[ij]-1),
                                              "_",substr(location,spaces[ij]+1,
                                              nchar(location)),sep="") else
          if (spaces[ij]==1) location <- substr(location,2,nchar(location))
      }

      apos <- instring("'",location)
      if (length(apos)>0) {
         for (ij in 1:length(apos))
           location <- paste(substr(location,1,apos[ij]-1),
                          substr(location,apos[ij]+1,nchar(location)),sep="")
      }

      slash.hunt <- instring("/",location)
      if (length(slash.hunt)>0) {
        for (iii in 1:length(slash.hunt)) {
          if (slash.hunt[iii]>0) location <-
            paste(substr(location,1,slash.hunt[1]-1),
            substr(location,slash.hunt[1]+1,nchar(location)),sep="-")
        }
      }
      name[i] <- location
  

  paren1 <- instring("(",header[2]); paren2 <- instring(")",header[2]); 
  country <- substr(header[2],paren1+1,paren2-1)
  i1n <- instring("N,",header[3])[1]
  i1s <- instring("S,",header[3])[1]
  i1 <- max(c(i1n,i1s),na.rm=TRUE)
  i2e <- instring("E,",header[3])[1]
  i2w <- instring("W,",header[3])[1]
  i2 <- max(c(i2e,i2w),na.rm=TRUE)
  i3 <-  instring("m ",header[3])[1]
  lat <- as.numeric(substr(header[3],15,i1-1))
  lon <- as.numeric(substr(header[3],i1+3,i2-1))
  alt <- as.numeric(substr(header[3],i2+3,i3-1))
  icode <- instring("code:",header[4])
  station <- as.numeric(strip(substr(header[4],icode+5,nchar(header[4]))))
  col.names <- c("YEAR","JAN","FEB","MAR","APR","MAY","JUN",
                   "JUL","AUG","SEP","OCT","NOV","DEC")
  ref <- paste("ClimateExplorer:",substr(header[5],3,nchar(header[5])))
  t2m <- read.table(fullURL,col.names=col.names)
  x <- cbind(t2m$JAN,t2m$FEB,t2m$MAR,t2m$APR,t2m$MAY,t2m$JUN,
             t2m$JUL,t2m$AUG,t2m$SEP,t2m$OCT,t2m$NOV,t2m$DEC)
  x[abs(x) >= 999] <- NA

  yy <- t2m$YEAR
  doubles <- yy[as.numeric(table(yy))>1]
  if (length(doubles)>0) {
    for (id in 1:length(doubles)) {
      iii <- (1:length(yy))[is.element(yy,doubles[id])]
      if (sum(is.finite(x[iii[2],])) > sum(is.finite(x[iii[1],]))) yy[iii][1] <- NA else
                                                                   yy[iii][2] <- NA
    }
    keep <- is.finite(yy)
    yy <- yy[keep]
    x <- x[keep,]
    srt <- order(yy)
    yy <- yy[srt]
    x <- x[srt,]
  }
  obs <- station.obj(x=x,yy=yy,
                     obs.name="Temperature",unit="deg C",ele=101,
                     station=station,lat=lat,lon=lon,alt=alt,
                     location=location,country=country,
                     ref=ref)
  invisible(obs)  
}
