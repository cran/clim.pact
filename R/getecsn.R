getecsn <- function(location="prompt",param=c("TG","RR"),dataset="blended",
                    data.path="data.eca",country=NULL,update=FALSE) {

  countries    <- c("Austria","Belgium","Denmark","Finland","France","Germany","Greece","Hungary",
                    "Iceland","Ireland","Italy","Luxembourg","Netherlands","Norway","Portugal",
                    "Spain","Sweden","Switzerland","United Kingdom","Russia","Greece",
                    "Romania","Israel","Algerie","Lituania","Ukraine","Yugoslavia","Estonia","Bosnia",
                    "Albania","Armenia","Greenland")
  country.code <- c("AT","BE","DK","FI","FR","DE","GR","HU",
                    "IS","IE","IT","LU","NL","NO","PT",
                    "ES","SE","CH","GB","RU","GR","IE",
                    "RO","IL","DZ","LT","UA","YU","EE","BA",
                    "AL","AM","GL")
  create.paths<-function(data.path) {
    if(!file.exists(data.path)) {
       if(!file.exists(dirname(data.path))) {
           create.paths( dirname(data.path) )
       }
    dir.create(data.path)
    }
  }

  print("getecsn")
  url.file <- switch(lower.case(dataset),
       "blended"="http://eca.knmi.nl/utils/downloadfile.php?file=download/ECA_blend_all.zip",
       "gcos"="http://eca.knmi.nl/utils/downloadfile.php?file=download/GCOS_blend_all.zip",
       "eumetnet"="http://eca.knmi.nl/utils/downloadfile.php?file=download/EUMETNET_blend_all.zip")
  #print(url.file)
  dash <- instring("/",url.file)
  f.name <- substr(url.file,dash[5]+1,nchar(url.file))
  #print(f.name)
  wd.0 <- getwd()
  print(wd.0)
  if (!file.exists(data.path)) create.paths(data.path)
  setwd(data.path)
  if ((!file.exists(f.name)) | update) {
    print(paste("download file from URL & store locally in",data.path))
    print("please wait a little while - download takes some time")
    download.file(url.file,f.name,method="internal")
    
    print("unzip file..."); print(getwd()); print(list.files(pattern=".zip")); 
    print(f.name)
    nc <- nchar(f.name)
#  zz <- unz(f.name,filename="location.txt",open="r") 
    zz <- unz(f.name,filename="location.txt") 
    system("unzip *.zip")
  
    a <- readLines(zz)
    close(zz)
    print("metadata OK")
    a <- a[20:length(a)]
    for (i in 1:length(a)) {
      while (length(instring(",",a[i])) < 5) { a[i] <- paste(a[i],",NA", sep="")}
      bad.char <- instring("^",a[i])
      if (length(bad.char)>0) {
      #print(a[i])
        a[i] <- paste(substr(a[i],1,bad.char-1),substr(a[i],bad.char+1,nchar(a[i])),sep="")
      #print(a[i])
      }
    }
    #print("save in 'eca.loc.txt'")
    writeLines(a,"eca.loc.txt")
  }
  #print("read as table")
  locs <- read.table("eca.loc.txt", sep =",",header=FALSE,
                     col.names=c("LOCID","LOCNAME","CN","LAT","LON","HGHT"))

  #print(summary(locs))
  n <- nchar(location)
  print(" get the data...")

  if (location!="prompt") {
    #print(location)
    loc <- lower.case(location)
    match <- is.element(substr(lower.case(as.character(locs$LOCNAME)),1,n),loc)
    #print("match")
  } else if (!is.null(country)) {
    print("country")
    contr <- lower.case(substr(countries,1,nchar(country)))
    ccode <- country.code[is.element(contr,lower.case(country))]
    match <- is.element(locs$CN,ccode)
    #print("match")
  }  else if ((location=="prompt") & (is.null(country))) match <- rep(TRUE,length(locs$LOCNAME))

  #print((1:length(locs$LOCNAME))[match])
  if (sum(match)==0) stop(paste('Could not find',location)) else
   if (sum(match)>1) {
    print(paste("Found",sum(match),"matches:")) 
    for (i in 1:sum(match)) print(paste(i,locs$LOCNAME[match][i]))
    ii <- as.numeric(readline('Type in the number of desired location: '))
   } else ii <- 1
  loc.id <- as.character(locs$LOCID[match][ii])
  station <- as.numeric(loc.id)
  lat.str <- as.character(locs$LAT[match][ii])
  lon.str <- as.character(locs$LON[match][ii])
  location <-  as.character(locs$LOCNAME[match][ii])
  col <- instring(":",lat.str)
  lat <- as.numeric(substr(lat.str,1,col[1]-1)) +
         100/60*as.numeric(substr(lat.str,col[1]+1,col[2]-1)) +
         100/60*100/60*as.numeric(substr(lat.str,col[2]+1,col[3]-1))
  col <- instring(":",lon.str)
  lon <- as.numeric(substr(lon.str,1,col[1]-1)) +
         100/60*as.numeric(substr(lon.str,col[1]+1,col[2]-1)) +
         100/60*100/60*as.numeric(substr(lon.str,col[2]+1,col[3]-1))
  alt <- locs$HGHT[match][ii]
  loc.id <- as.character(as.numeric(loc.id))
  #print(loc.id)
  while (nchar(loc.id) < 6) loc.id <- paste("0",loc.id,sep="")
  #print(loc.id)
  loc.name <- paste("LOCID",loc.id,sep="")
  for (iele in 1:length(param)) {
    station.file <- paste(upper.case(param[iele]),"_",loc.name,".txt",sep="")
    print(paste("Read from",f.name,"file called",station.file))
    #in.data <- unz(f.name,filename=station.file,open="r")
    #a <- readLines(in.data)
    #close(in.data)
    in.data <- file(station.file,open="r")
    a <- readLines(in.data)
    #print(length(a))
    ref <- paste(a[3],a[4],a[5])
    #print(ref)
    country <- a[15]
    #print(country)
    a <- a[20:length(a)]
    #print("save in 'eca.data.txt'")
    writeLines(a[1:length(a)],"eca.data.txt")
#    eval(parse(text=paste("data.",i,sep=""))) <- read.table("eca.data.txt", sep =",",header=FALSE)    
    if (iele==1) data.1 <- read.table("eca.data.txt", sep =",",header=FALSE)    
    if (iele==2) data.2 <- read.table("eca.data.txt", sep =",",header=FALSE)    
    #print(summary(eval(parse(text=paste("data.",i,sep="")))))
  }
  
  #print(lower.case(country))
  for (i in 1:length(countries)) if (length(grep(lower.case(countries[i]),lower.case(country)))>=1) country <- countries[i]
  #print(country)

  i1 <- is.element(data.1$V2,data.2$V2)
  i2 <- is.element(data.2$V2,data.1$V2)
 
  yy1 <- floor(data.1$V2/10000); mm1 <- floor((data.1$V2 - yy1*10000)/100); dd1 <- data.1$V2 - yy1*10000 - mm1*100
  yy2 <- floor(data.2$V2/10000); mm2 <- floor((data.2$V2 - yy2*10000)/100); dd2 <- data.2$V2 - yy2*20000 - mm2*100

  yy <- yy1[i1]; mm <- mm1[i1]; dd <- dd1[i1]
  t2m <- data.1$V3[i1]/10;    t2m[t2m <= -99] <- NA
  precip <- data.2$V3[i2]/10; precip[precip <= -99] <- NA
  setwd(wd.0)

  obs.name <- c(" "," "); unit <- obs.name
  obs.name[1] <- switch(lower.case(param[1]),"tg"="Daily mean temperature","tx"="Daily max temperature",
                                   "tn"="daily minumum temperature","rr"="24-hour precipitation",
                                   "pp"="daily mean sea level pressure")
  obs.name[2] <- switch(lower.case(param[2]),"tg"="Daily mean temperature","tx"="Daily max temperature",
                                   "tn"="daily minumum temperature","rr"="24-hour precipitation",
                                   "pp"="daily mean sea level pressure")
  unit[1] <- switch(lower.case(param[1]),"tg"="deg C","tx"="deg C",
                                   "tn"="deg C","rr"="mm","pp"="hPa")
  unit[2] <- switch(lower.case(param[2]),"tg"="deg C","tx"="deg C",
                                   "tn"="deg C","rr"="mm","pp"="hPa")

  ecsn <- station.obj.dm(t2m=t2m,precip=precip,dd=dd,mm=mm,yy=yy,
                         obs.name=obs.name,unit=unit,ele=param,
                         station=station,lat=lat,lon=lon,alt=alt,
                         location=location,
                         country=country,ref=ref)
  invisible(ecsn)
}
