# R.E. Benestad, met.no, Oslo, Norway 16.04.2002
# rasmus.benestad@met.no
#
# Modified 27.04.2004 to base the IO on the ncdf package instead of
# netCDF (which is being phased out).
#------------------------------------------------------------------------


retrieve.nc <- function(filename="data/ncep_t2m.nc",v.nam="AUTO",
                        l.scale=TRUE,greenwich=TRUE,
                        x.nam="lon",y.nam="lat",z.nam="lev",t.nam="tim",
                        x.rng=NULL,y.rng=NULL,t.rng=NULL,force.chron=TRUE) {
  library(ncdf)
  if (!file.exists(filename)) {
    stop(paste("Sorry,",filename," does not exist!"))
  }

  dat.att <- cdfcont(filename)
  ncid1 <- open.ncdf(filename)
  if (v.nam=="AUTO") {
    v1 <- ncid1$var
    n.vars <- length(v1)
    if (n.vars==1) v1 <- v1[[1]] else {
      ipick<-0
      for (i in seq(n.vars,1,by=-1)) {
        if (v1[[i]]$ndim==3) ipick <- i
      }
      if (ipick > 0) v1 <- v1[[ipick]] else {
         print("Tip: use open.ncdf to read the data manually")
         print(names(v1))
         print(dat.att)
         stop("Error: did'n find a variable with 3D")
      }
    }
  } else {
    v1 <- ncid1$var
    vars <- names(v1)
    ipick <- grep(v.nam,vars)
    if (length(ipick)==0) {
      print(vars)
      ipick <- as.numeric(readline(paste("Choose variable (1 - ",length(vars),"): ",sep="")))
    }
    v1 <- ncid1$var[[ipick]]
  }
  data <- get.var.ncdf(ncid1,v1)
  close.ncdf(ncid1)
  vars <- v1$name
  nvars <-  ncid1$nvars
  print(paste("Reading",vars))
  v.nam <- v1$longname
  dims <- names(ncid1$dim)
  vars <- names(ncid1$dim)
  n.dim <- ncid1$ndims
  d <- rep(0,nvars)
  dat <- NULL

  ilon <- grep("lon",lower.case(names(ncid1$dim)))
  ilat <- grep("lat",lower.case(names(ncid1$dim)))
  itim <- grep("tim",lower.case(names(ncid1$dim)))
  ilev <- grep("lev",lower.case(names(ncid1$dim)))
  #print(c(ilon,ilat,itim))

  eval(parse(text=paste("lon <- ncid1$dim$",names(ncid1$dim)[ilon],"$vals",sep="")))
  eval(parse(text=paste("lat <- ncid1$dim$",names(ncid1$dim)[ilat],"$vals",sep="")))
  eval(parse(text=paste("tim <- ncid1$dim$",names(ncid1$dim)[itim],"$vals",sep="")))
  attr(lon,"unit") <- eval(parse(text=paste("ncid1$dim$",names(ncid1$dim)[ilon],"$units",sep="")))
  attr(lat,"unit") <- eval(parse(text=paste("ncid1$dim$",names(ncid1$dim)[ilat],"$units",sep="")))
  attr(tim,"time_origin") <- dat.att$time.origin
  if (!is.null(dat.att$time.unit)) attr(tim,"unit") <- dat.att$time.unit else 
    attr(tim,"unit") <- eval(parse(text=paste("ncid1$dim$",names(ncid1$dim)[itim],"$units",sep="")))
   
  if (length(ilev)>0) {
    eval(parse(text=paste("lev <- ncid1$dim$",names(ncid1$dim)[ilev],"$vals",sep="")))
    attr(lev,"unit") <- eval(parse(text=paste("ncid1$dim$",names(ncid1$dim)[ilev],"$units",sep="")))
    n.dim <- 4
  } else {
    lev <- NULL
    n.dim <- 3
  }
 # Re-order the data: (old convention)
  nt <- length(tim); ny <- length(lat); nx <- length(lon)
  #print(c(nt,ny,nx,NA,dim(data)))
  dat <- data*NA; dim(dat) <- c(nt,ny,nx)
  for (i in 1:nt) dat[i,,] <- t(as.matrix(data[,,i]))
  
  cmon<-c('Jan','Feb','Mar','Apr','May','Jun',
          'Jul','Aug','Sep','Oct','Nov','Dec')
  season<-cbind(c(12,1,2),c(3,4,5),c(6,7,8),c(9,10,11))
  season.c<-c("","DJF","MAM","JJA","SON")

  dtim <- diff(tim)
  if ( sum(dtim<=0) > 0) {
    print(paste("Warning! Test of chonological order finds",sum(dtim<=0),"jump(s)"))
    print(paste("median(dtim)=",median(dtim)))
    if (force.chron) {
      nt <- length(tim)
      tim.att <- attributes(tim)
      dtims <- as.numeric(row.names(table(dtim)))
      if (length(dtims < 4)) {
        print(paste("Force correction: assume tim[1] is correct,",
                    median(dtim),"is correct time step, and length=",nt))
        tim <- seq(tim[1],tim[1]+nt-1,by=median(dtim))
      } else {
        dt <- readline("What is the correct time step? (0 leaves tim unchanged)")
        if (dt != 0) tim <- seq(tim[1],tim[1]+nt-1,by=dt)
      }
    }
    print(paste("length(tim)=",length(tim),"nt=",nt))
#    print("set new attributes for tim")
    attributes(tim) <- tim.att
#    print("continue...")
  }

  t.unit <- attr(tim,"unit")
  dat.att$unit <-v1$units

  if (!is.null(dat.att$time.origin)) {
    torg <-  dat.att$time.origin
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
  if (yy0==0) {
    print('There is no year zero (Press et al., Numerical recipies)')
    print("'> print(julday(1,1,1)-julday(1,1,-1))' gives 365")
    print('julday wont work unless the time is fixed')
    print("year0 is set to 1, and 365 days is subtracted from tim")
    if (substr(lower.case(t.unit),1,4)=="hour") tim <- tim - 365*24
    if (substr(lower.case(t.unit),1,3)=="day") tim <- tim - 365
    if (substr(lower.case(t.unit),1,3)=="mon") tim <- tim - 12
    if (substr(lower.case(t.unit),1,5)=="year") tim <- tim - 1
    yy0 <- 1
  }

  print(paste("Time unit:",lower.case(t.unit)))
  if (substr(lower.case(t.unit),1,3)=="mon") {
    tim <- floor(tim)
    mm <- mod(mm0 + tim - 1,12)+1
    yy  <- yy0 + floor((tim+mm0-1)/12)
    dd <- rep(15,length(tim))
    obj.type <- "monthly.field.object"
  } else if (substr(lower.case(t.unit),1,3)=="day") {
    mmddyy <- caldat(tim + julday(mm0,dd0,yy0))
    mm <- mmddyy$month
    yy <- mmddyy$year
    dd <- mmddyy$day
    obj.type <- "daily.field.object"
  } else if (substr(lower.case(t.unit),1,4)=="hour") {
    mmddyy <- caldat(tim/24 + julday(mm0,dd0,yy0))
    mm <- mmddyy$month
    yy <- mmddyy$year
    dd <- mmddyy$day
    t.unit <- "day"
    obj.type <- "field.object"
  } 
  
# Extra processing for NCEP files e.g. with Time unit: hours since 1-1-1 00:00:0.0.
  if ( ((substr(lower.case(t.unit),1,4)=="hour") |
        (substr(lower.case(t.unit),1,3)=="day")) &
       (max(diff(dd)) == 0) ) {
    print("Monthly data, but time unit set to 'hour'/'day'")
    print("Set time unit to month")
    obj.type <- "monthly.field.object"
    dd[] <- 15
    t.unit <- "month"
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
    torg <- attr(tim,"time.origin")
    tunit <- attr(tim,"unit")
    tim <- tim[t.keep]
    attr(tim,"time.origin") <- torg
    attr(tim,"unit") <- tunit
    yy <- yy[t.keep]
    mm <- mm[t.keep]
    dd <- dd[t.keep]
    nt <- length(tim)
  }
  print(paste("First & last records:",yy[1],mm[1],dd[1],
              "&",yy[length(yy)],mm[length(mm)],dd[length(dd)]))
  
#  print(dat.att$scale.factor)
#  print(dat.att$add.offset)
#  print(dat.att$unit)

#  print(dat.att)
  if ((l.scale) & !is.null(dat.att$scale.factor)) {
     if (is.finite(dat.att$scale.factor)) dat <- dat * dat.att$scale.factor
  }
  # Have included a sanity test to detect an old 'bug': offset 273 and
  # units of deg C..
  if ( ((l.scale) & !is.null(dat.att$add.offset))) {
      if ( (dat.att$add.offset!=273) &
           (dat.att$unit=="deg C")) {
        a <- readline(prompt="Correct an old bug? (y/n)")
        if (lower.case(a)=="y") dat <- dat + dat.att$add.offset} else
        if (is.finite(dat.att$add.offset)) dat <- dat + dat.att$add.offset
  }
  if (l.scale) {   
    print("BEFORE scale adjustment & weeding")
    print(summary(as.vector(dat)))
    dat[dat == dat.att$missing.value] <- NA
    if (sum(is.na(dat))>0) print(paste(sum(is.na(dat)),"of",length(dat),
                                 " are set to 'NA'"))
      print("AFTER scale adjustment & weeding")
  }

  if (!is.null(dat.att$units)) {
     if (is.finite(dat.att$units)) dat.att$unit <- dat.att$units
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
  slash <- instring("/",filename)
  dot <- instring(".",filename)

  id.t <- rep(substr(filename,slash[length(slash)]+1,
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

