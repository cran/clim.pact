retrieve.nc <- function(f.name="data/ncep_t2m.nc",v.nam="AUTO",
                        l.scale=TRUE,greenwich=TRUE,
                        x.nam="lon",y.nam="lat",t.nam="tim",
                        x.rng=NULL,y.rng=NULL,t.rng=NULL) {
  library(netCDF)
  library(chron)
  if (!file.exists(f.name)) {
    stop(paste("Sorry,",f.name," does not exist!"))
  }
  ncid1<-open.netCDF(f.name)
  dat <- read.netCDF(ncid1)
  close.netCDF(ncid1)
  vars <- names(dat)
  nvars <- length(vars)
  d <- rep(0,nvars)
#  print("Searching for variables")
  for (i in 1:nvars) {
#      print(vars[i])
      if (sum(grep("_",vars[i]))==0) {
      expr <- parse(text=paste("dat$",vars[i],sep=""))
      if (lower.case(substr(vars[i],1,nchar(x.nam)))==x.nam) lon <- eval(expr)
      expr <- parse(text=paste("dat$",vars[i],sep=""))
      if (lower.case(substr(vars[i],1,nchar(y.nam)))==y.nam) lat <- eval(expr)
      expr <- parse(text=paste("dat$",vars[i],sep=""))
      if (lower.case(substr(vars[i],1,nchar(t.nam)))==t.nam) tim <- eval(expr)

      expr <- parse(text=paste("dim(dat$",vars[i],")",sep=""))
      size <- eval(expr)
      expr <- parse(text=paste("dat$",vars[i],sep=""))
      if (lower.case(substr(vars[i],1,nchar(v.nam)))==t.nam) {
        dat <- eval(expr)
      } else if ((length(size)==3) & (v.nam=="AUTO")) {
        dat <- eval(expr)
        var.name <- vars[i]
      } else if (v.nam=="ASK") {
        for (ii in 1:nvars) {
          print(vars[i])
          i.var <- readline(prompt="Select variable:")
          expr <- parse(text=paste("dat$",vars[i.var],sep=""))
          dat <- eval(expr)
          var.name <- vars[i.var]
        }
      }
    }
  }
  slash <- instring("/",f.name)
  dot <- instring(".",f.name)
  nx <- length(lon)
  ny <- length(lat)
  nt <- length(tim)
  dat.att <- attributes(dat)
  
  cmon<-c('Jan','Feb','Mar','Apr','May','Jun',
          'Jul','Aug','Sep','Oct','Nov','Dec')
  season<-cbind(c(12,1,2),c(3,4,5),c(6,7,8),c(9,10,11))
  season.c<-c("","DJF","MAM","JJA","SON")

  t.unit <- NULL
  if (!is.null(attributes(tim)$unit)) t.unit <- attributes(tim)$unit else
    if (!is.null(attributes(tim)$units)) t.unit <- attributes(tim)$units
  
  print("Time information:")
  torg <-  attributes(tim)$"time_origin"
  if (!is.null(torg)) {
    yy0 <- as.numeric(substr(torg,8,11))
    dd0 <- as.numeric(substr(torg,1,2))
    mm0 <- switch(substr(torg,4,6),
                  "Jan"=1,"Feb"=2,"Mar"=3,"Apr"=4,"May"=5,"Jun"=6,
                  "Jul"=7,"Aug"=8,"Sep"=9,"Oct"=10,"Nov"=11,"Dec"=12)
  } else if (grep("since",lower.case(t.unit))) {
    # Format: time:units = "hours since 1-1-1 00:00:0.0" (NCEP reanalysis)
    t.org.pos <- regexpr("since",lower.case(t.unit))
    torg  <- substr(t.unit,t.org.pos+6,nchar(t.unit))
    dash <- instring("-",torg)
    spc <- instring(" ",torg)
    yy0 <- as.numeric(substr(torg,1,dash[1]-1))
    mm0 <- as.numeric(substr(torg,dash[1]+1,dash[2]-1))
    dd0 <- as.numeric(substr(torg,dash[2]+1,spc[1]-1))
  }
  print(c(yy0,mm0,dd0))
  
  print(paste("Time unit:",lower.case(t.unit)))
  if (substr(lower.case(t.unit),1,5)=="month") {
    tim <- floor(tim)
    mm <- mod(mm0 + tim - 1,12)+1
    yy  <- yy0 + floor((tim+mm0-1)/12)
    dd <- rep(15,length(tim))
    obj.type <- "monthly.field.object"
  } else if (substr(lower.case(t.unit),1,3)=="day") {
    mmddyy<-month.day.year(tim,origin=c(mm0,dd0,yy0))
    mm <- mmddyy$month
    yy <- mmddyy$year
    dd <- mmddyy$day
    obj.type <- "daily.field.object"
  } else if (substr(lower.case(t.unit),1,4)=="hour") {
    mmddyy<-month.day.year(tim/24,origin=c(mm0,dd0,yy0))
    mm <- mmddyy$month
    yy <- mmddyy$year
    dd <- mmddyy$day
    obj.type <- "field.object"
  }
#  print("Latitude:")
  if (attributes(lat)$"unit"=="degrees_south") lat <- lat * -1
  if (attributes(lon)$"unit"=="degrees_west") lon <- lon * -1
  if (!is.null(y.rng)) {
    print(range(lat))
    print("Extract latitudes:")
    print(y.rng)
    y.keep <- (lat >= min(y.rng)) & (lat <= max(y.rng))
    dat <- dat[,y.keep,]
    lat <- lat[y.keep]
  }
  if (greenwich) {
    lon[lon > 180] <- lon[lon > 180]-360
  }
#  print("Sort longs and lats")
  x.srt <- order(lon)
  y.srt <- order(lat)
  lon <- lon[x.srt]
  lat <- lat[y.srt]
  dat <- dat[,y.srt,x.srt]
  if (!is.null(x.rng)) {
    print(range(lon))
    print("Extract longitudes:")
    print(x.rng)
    x.keep <- (lon >= min(x.rng)) & (lon <= max(x.rng))
    dat <- dat[,,x.keep]
    lon <- lon[x.keep]
  }
  if (!is.null(t.rng)) {
    print(range(yy))
    print("Extract times:")
    print(t.rng)
    t.keep <- (yy >= min(t.rng)) & (yy <= max(t.rng))
    dat <- dat[t.keep,,]
    torg <- attr(tim,"time_origin")
    tunit <- attr(tim,"unit")
    tim <- tim[t.keep]
    attr(tim,"time_origin") <- torg
    attr(tim,"unit") <- tunit
    yy <- yy[t.keep]
    mm <- mm[t.keep]
    dd <- dd[t.keep]
    nt <- length(tim)
  }
  print(paste("First & last records:",yy[1],mm[1],dd[1],
              "&",yy[length(yy)],mm[length(mm)],dd[length(dd)]))
#  print(attributes(dat)$"scale_factor")
#  print(attributes(dat)$"add_offset")
#  print(attributes(dat)$"unit")

#  print(dat.att)
  if (l.scale) {
    dat <- dat * dat.att$"scale_factor"
  }
  # Have included a sanity test to detect an old 'bug': offset 273 and
  # units of deg C..
  if ( (l.scale) |
      (dat.att$"add_offset"!=273) &
      (dat.att$"unit"=="deg C")) {
     dat <- dat + dat.att$"add_offset"}
  if (l.scale) {   
    print("BEFORE scale adjustment & weeding")
    print(summary(as.vector(dat)))
    dat[dat == dat.att$"missing_value"] <- NA
    if (sum(is.na(dat))>0) print(paste(sum(is.na(dat)),"of",length(dat),
                                 " are set to 'NA'"))
      print("AFTER scale adjustment & weeding")
  }

  if ((dat.att$unit=="K") | (dat.att$unit=="Kelvin") |
      (dat.att$unit=="degrees Kelvin") |
      (dat.att$unit=="deg K")) {
    dat <- dat - 273
    dat.att$unit <- "deg C"
  }
    if ((dat.att$unit=="Pa") | (dat.att$unit=="Pascal") |
      (dat.att$unit=="N/m^2") |
      (dat.att$unit=="N m^{-1}")) {
    dat <- dat/100
    dat.att$unit <- "hPa"
  }
  print(summary(as.vector(dat)))
  nx <- length(lon)
  ny <- length(lat)
  print(c(nt,ny,nx))
  eos <- nchar(var.name)
  if (instring("-",var.name)> 0) {
    eos <- instring("-",var.name)-1
  } else if (instring("_",var.name)> 0) {
    eos <- instring("_",var.name)-1
  }
  var.name <- substr(var.name,1,eos)
  id.x <- matrix(rep(var.name,ny*nx),ny,nx)
  id.t <- rep(substr(f.name,slash[length(slash)]+1,
                     dot[length(dot)]-1),nt)              
  
  retrieve.nc  <- list(dat=dat,lon=lon,lat=lat,tim=tim,v.name=var.name,
                       id.x=id.x,id.t=id.t,yy=yy,mm=mm,dd=dd,n.fld=1,
                       id.lon=rep(var.name,nx),id.lat=rep(var.name,ny))
  class(retrieve.nc) <- c("field",obj.type)
  invisible(retrieve.nc)
}

