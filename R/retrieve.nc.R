# R.E. Benestad, met.no, Oslo, Norway 16.04.2002
# rasmus.benestad@met.no
#------------------------------------------------------------------------


retrieve.nc <- function(f.name="data/ncep_t2m.nc",v.nam="AUTO",
                        l.scale=TRUE,greenwich=TRUE,
                        x.nam="lon",y.nam="lat",z.nam="lev",t.nam="tim",
                        x.rng=NULL,y.rng=NULL,t.rng=NULL) {
  library(netCDF)
#  library(chron)
#  library(date)
  if (!file.exists(f.name)) {
    stop(paste("Sorry,",f.name," does not exist!"))
  }
  ncid1<-open.netCDF(f.name)
  data <- read.netCDF(ncid1)
  close.netCDF(ncid1)
  vars <- names(data)
  nvars <- length(vars)
  n.dim <- 3
  if (length(grep(z.nam,lower.case(substr(vars,1,3))))==1) {
    n.dim <- 4
  } else {
    lev <- NULL
    nz <- 0
  }
  d <- rep(0,nvars)
  dat <- NULL
  
#  print("Searching for variables")
  for (i in 1:nvars) {
      expr <- parse(text=paste("dim(data$'",vars[i],"')",sep=""))
      size <- eval(expr)

      expr <- parse(text=paste("lon <- data$'",vars[i],"'",sep=""))
      if (lower.case(substr(vars[i],1,nchar(x.nam)))==lower.case(x.nam) &
          (length(size)==1)) {
        eval(expr)
      }
      expr <- parse(text=paste("lat <- data$'",vars[i],"'",sep=""))
      if (lower.case(substr(vars[i],1,nchar(y.nam)))==lower.case(y.nam) &
          (length(size)==1)) {
        eval(expr)
      }
      expr <- parse(text=paste("tim <- data$'",vars[i],"'",sep=""))
      if (lower.case(substr(vars[i],1,nchar(t.nam)))==lower.case(t.nam) &
          (length(size)==1)) {
        eval(expr)
      }

      if (n.dim==4) {
        expr <- parse(text=paste("lev <- data$'",vars[i],"'",sep=""))
        if (lower.case(substr(vars[i],1,nchar(z.nam)))==lower.case(z.nam) &
          (length(size)==1)) {
          eval(expr)
        }
      }
#      print(paste('dim(data$',vars[i],')',sep=""))
#      print(size)
      expr <- parse(text=paste("dat <- data$'",vars[i],"'",sep=""))
      if ((lower.case(substr(vars[i],1,nchar(v.nam)))==lower.case(v.nam)) &
          (is.null(dat))) {
        eval(expr)
        if (n.dim==3) print(paste("Data:",vars[i]," dim:",
                       dim(dat)[1],"x",dim(dat)[2],"x",dim(dat)[3])) else
                      print(paste("Data:",vars[i]," dim:",
                         dim(dat)[1],"x",dim(dat)[2],"x",
                                  dim(dat)[3],"x",dim(dat)[4]))
      } else if ((length(size)==n.dim) & (v.nam=="AUTO")) {
        eval(expr)
        v.nam <- vars[i]
        print(paste("Data:",v.nam," dim:",
                    dim(dat)[1],"x",dim(dat)[2],"x",dim(dat)[3]))
      } else if (v.nam=="ASK") {
          print(vars)
          i.var <- as.numeric(readline(prompt=
                   paste("Select data variable (1-",nvars,"): ",sep="")))
          print(paste("dat <- data$",vars[i.var],sep=""))
          expr <- parse(text=paste("dat <- data$'",vars[i.var],"'",sep=""))
          eval(expr)
          v.nam <- vars[i.var]
          print(paste("Data:",v.nam," dim:",
                dim(dat)[1],"x",dim(dat)[2],"x",dim(dat)[3]))
      } 
  }

  if (is.null(dat)) {
    print("Did not find the data")
    print(vars)
    print("Try the option v.nam='ASK'")
    return()
  }
#  print("Found all variables")
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
  if (!is.null(attributes(tim)$"time_origin")) {
    torg <-  attributes(tim)$"time_origin"
  } else torg <- NULL

  if (!is.null(torg)) {
    yy0 <- as.numeric(substr(torg,8,11))
    dd0 <- as.numeric(substr(torg,1,2))
    mm0 <- switch(lower.case(substr(torg,4,6)),
                  "jan"=1,"feb"=2,"mar"=3,"apr"=4,"may"=5,"jun"=6,
                  "jul"=7,"aug"=8,"sep"=9,"oct"=10,"nov"=11,"dec"=12)
  } else if (grep("since",lower.case(t.unit))) {
    # Format: time:units = "hours since 1-1-1 00:00:0.0" (NCEP reanalysis)
    t.org.pos <- regexpr("since",lower.case(t.unit))
    torg  <- substr(t.unit,t.org.pos+6,nchar(t.unit))
    print(paste("torg=",torg))
    dash <- instring("-",torg)
    spc <- instring(" ",torg)
    yy0 <- as.numeric(substr(torg,1,dash[1]-1))
    mm0 <- as.numeric(substr(torg,dash[1]+1,dash[2]-1))
    dd0 <- as.numeric(substr(torg,dash[2]+1,spc[1]-1))
    if (is.na(dd0)) dd0  <- 15
  }
  print(paste("Time origin: (year-month-day)",yy0,"-",mm0,"-",dd0))
  
  print(paste("Time unit:",lower.case(t.unit)))
  if (substr(lower.case(t.unit),1,5)=="month") {
    tim <- floor(tim)
    mm <- mod(mm0 + tim - 1,12)+1
    yy  <- yy0 + floor((tim+mm0-1)/12)
    dd <- rep(15,length(tim))
    obj.type <- "monthly.field.object"
  } else if (substr(lower.case(t.unit),1,3)=="day") {
#    mmddyy<-month.day.year(tim,origin=c(mm0,dd0,yy0))
#    mmddyy <- date.mdy(tim + mdy.date(mm0,dd0,yy0,nineteen=FALSE), weekday = FALSE)
    mmddyy <- caldat(tim + julday(mm0,dd0,yy0))
    mm <- mmddyy$month
    yy <- mmddyy$year
    dd <- mmddyy$day
    obj.type <- "daily.field.object"
  } else if (substr(lower.case(t.unit),1,4)=="hour") {
#    mmddyy<-month.day.year(tim/24,origin=c(mm0,dd0,yy0))
#    mmddyy <- date.mdy(tim/24 + mdy.date(mm0,dd0,yy0,nineteen=FALSE),weekday = FALSE)
    mmddyy <- caldat(tim/24 + julday(mm0,dd0,yy0))
    mm <- mmddyy$month
    yy <- mmddyy$year
    dd <- mmddyy$day
    t.unit <- "day"
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
    if (n.dim==3) dat <- dat[,y.keep,] else
                  dat <- dat[,,y.keep,] 
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
  if (n.dim==3) dat <- dat[,y.srt,x.srt] else
                dat <- dat[,,y.srt,x.srt]
  
  if (!is.null(x.rng)) {
    print(range(lon))
    print("Extract longitudes:")
    print(x.rng)
    x.keep <- (lon >= min(x.rng)) & (lon <= max(x.rng))
    if (n.dim==3) dat <- dat[,,x.keep] else
                  dat <- dat[,,,x.keep]
    lon <- lon[x.keep]
  }
  if (!is.null(t.rng)) {
    print(range(yy))
    print("Extract times:")
    print(t.rng)
    t.keep <- (yy >= min(t.rng)) & (yy <= max(t.rng))
    if (n.dim==3) dat <- dat[t.keep,,] else
                  dat <- dat[t.keep,,,]
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
  if ((l.scale) & !is.null(attributes(dat)$"scale_factor")) {
    dat <- dat * dat.att$"scale_factor"
  }
  # Have included a sanity test to detect an old 'bug': offset 273 and
  # units of deg C..
  if ( ((l.scale) & !is.null(attributes(dat)$"add_offset"))) {
      if ( (dat.att$"add_offset"!=273) &
           (dat.att$"unit"=="deg C")) {
        a <- readline(prompt="Correct an old bug? (y/n)")
        if (lower.case(a)=="y") dat <- dat + dat.att$"add_offset"} else
        dat <- dat + dat.att$"add_offset"
  }
  if (l.scale) {   
    print("BEFORE scale adjustment & weeding")
    print(summary(as.vector(dat)))
    dat[dat == dat.att$"missing_value"] <- NA
    if (sum(is.na(dat))>0) print(paste(sum(is.na(dat)),"of",length(dat),
                                 " are set to 'NA'"))
      print("AFTER scale adjustment & weeding")
  }

  if (!is.null(dat.att$units)) {
     dat.att$unit <- dat.att$units
  } 
  if ((dat.att$unit=="K") | (dat.att$unit=="Kelvin") |
      (dat.att$unit=="degrees Kelvin") |
      (dat.att$unit=="deg K") | (dat.att$unit=="degK")) {
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
  eos <- nchar(v.nam)
  if (instring("-",v.nam)> 0) {
    eos <- instring("-",v.nam)-1
  } else if (instring("_",v.nam)> 0) {
    eos <- instring("_",v.nam)-1
  }
  v.nam <- substr(v.nam,1,eos)
  id.x <- matrix(rep(v.nam,ny*nx),ny,nx)
  id.t <- rep(substr(f.name,slash[length(slash)]+1,
                     dot[length(dot)]-1),nt)              
  dat.att$time.unit <- t.unit
  dat.att$time.origin <- torg
  retrieve.nc  <- list(dat=dat,lon=lon,lat=lat,tim=tim,lev=lev,
                       v.name=v.nam,id.x=id.x,id.t=id.t,
                       yy=yy,mm=mm,dd=dd,n.fld=1,
                       id.lon=rep(v.nam,nx),id.lat=rep(v.nam,ny),
                       attributes=dat.att)
  class(retrieve.nc) <- c("field",obj.type)
  invisible(retrieve.nc)
}

