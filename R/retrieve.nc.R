# R.E. Benestad, met.no, Oslo, Norway 16.04.2002
# rasmus.benestad@met.no
#
# Modified 27.04.2004 to base the IO on the ncdf package instead of
# netCDF (which is being phased out).
#------------------------------------------------------------------------


retrieve.nc <- function(filename=file.path("data","air.mon.mean.nc"),v.nam="AUTO",
                        l.scale=FALSE,greenwich=TRUE,silent=FALSE,
                        x.nam="lon",y.nam="lat",z.nam="lev",t.nam="tim",
                        x.rng=NULL,y.rng=NULL,z.rng=NULL,t.rng=NULL,
                        force.chron=TRUE,force365.25=FALSE,regular=TRUE,daysayear=365.25,
                        forceBC=TRUE,use.cdfcont=FALSE,torg=NULL,t.unit=NULL) {
#print("retrieve.nc:")
  library(ncdf)
  if (!file.exists(filename)) {
    stop(paste("Sorry,",filename," does not exist!"))
  }

  a <- Sys.info()
  ncid <- open.ncdf(filename)
  nv <- ncid$nvars
  cdfvars <- rep("-",nv)
  for (i in 1:nv) cdfvars[i] <- ncid$var[[i]]$name
  miss <-  ncid$var[[1]]$missval
 
  ipick <- 1
  if (nv > 1) {
    ipick <- grep(v.nam,cdfvars)
    if (length(ipick)==0) {
      print(cdfvars)
      ipick <- as.numeric(readline(paste("Choose variable (1 - ",length(cdfvars),"): ",sep="")))
    }
  } 
  v1 <- cdfvars[ipick]
 
  #nd <- ncid$ndims
  nd <- ncid$var[[ipick]]$ndims
  cdfdims <- rep("-",nd)
  for (i in 1:nd) cdfdims[i] <- ncid$var[[ipick]]$dim[[i]]$name
#print(cdfdims)
  ilon <- grep(x.nam,lower.case(cdfdims))
  ilat <- grep(y.nam,lower.case(cdfdims))
  itim <- grep(t.nam,lower.case(cdfdims))
  ilev <- grep(z.nam,lower.case(cdfdims))
#print(c(ilon,ilat,itim,ilev)) 
  scal <- NULL; offs <- NULL

  arv <- att.get.ncdf(ncid, cdfvars[ipick], 'scale_factor')
  if( arv$hasatt ) scal <- arv$value else
       scal <- 1

  arv <- att.get.ncdf(ncid, cdfvars[ipick], 'add_offset')
  if( arv$hasatt ) offs <- arv$value else
       offs <- 0

  arv <- att.get.ncdf(ncid, cdfvars[ipick], 'units')
  if( arv$hasatt ) unit <- arv$value else {
    arv <- att.get.ncdf(ncid, cdfvars[ipick], 'unit')
    if( arv$hasatt ) unit <- arv$value else
       print(paste("Attribute unit not found for",v1))  
       unit <- "unknown"
  }

  arv <- att.get.ncdf(ncid, cdfdims[itim], 'calendar')
  if( arv$hasatt ) calendar <- arv$value else
       calendar <- "ordinary"
print(calendar)
  if (calendar=="noleap") {
    print("Detected 'noleap' Calendar: set daysayear=365")
    daysayear <- 365
  }

  if ( (is.null(torg)) | is.null(t.unit)) {
    if ( (lower.case(a[1])=="linux") & (use.cdfcont) ) {
       if (!silent) print("Linux & use.cdfcont : call cdfconf()")
       if (is.null(torg)) torg <- cdfcont(filename)$time.origin
       if (is.null(t.unit)) t.unit <- cdfcont(filename)$time.unit
    } else {
       arv <- att.get.ncdf(ncid, cdfdims[itim], 'time_origin')
       if( arv$hasatt & is.null(torg)) torg <- arv$value else {
         if (is.null(torg)) print("Attribute time_origin not found")  
      }
     arv <- att.get.ncdf(ncid, cdfdims[itim], 'time_unit')
       if( arv$hasatt & is.null(t.unit)) t.unit <- arv$value else {
         arv <- att.get.ncdf(ncid, cdfdims[itim], 'unit')
         if( arv$hasatt & is.null(t.unit)) t.unit <- arv$value else {
           arv <- att.get.ncdf(ncid, cdfdims[itim], 'units')
           if( arv$hasatt & is.null(t.unit)) t.unit <- arv$value 
         } 
       }
    }
  }

#print(torg); print(t.unit)
  
  if (!silent) print(paste("Reading",v1))
  arv <- att.get.ncdf(ncid, v1, 'long_name')
  if( arv$hasatt ) lon.nam <- arv$value else lon.nam <- v1
  v.nam <- v1
  
  lon <- get.var.ncdf(ncid,cdfdims[ilon])
  lat <- get.var.ncdf(ncid,cdfdims[ilat])
  tim <- get.var.ncdf(ncid,cdfdims[itim])

#print("Attributes:")
  attr(lon,"unit") <- eval(parse(text=paste("ncid$dim$'",cdfdims[ilon],"'$units",sep="")))
  attr(lat,"unit") <- eval(parse(text=paste("ncid$dim$'",cdfdims[ilat],"'$units",sep="")))
  print(paste("ncid$dim$'",cdfdims[itim],"$units'",sep=""))
  attr(tim,"unit") <- eval(parse(text=paste("ncid$dim$'",cdfdims[itim],"'$units",sep="")))
  if (is.null(t.unit)) t.unit <- attr(tim,"unit")
  print(paste("Time, units: ",t.unit))

  if (length(ilev)>0) {
    print(paste("Get the levels: ",min(z.rng),max(z.rng)))
    lev <- get.var.ncdf(ncid,cdfdims[ilev])
    attr(lev,"unit") <- eval(parse(text=paste("ncid$dim$",cdfdims[ilev],"$units",sep="")))
# HERE 05.02.2008   
#    if (length(lev)==1) nd <- nd - 1
  } else {
    lev <- NULL
  }

  start <- rep(1,nd); count <- rep(1,nd)
# count <- c(length(lon),length(lat),length(tim))
  for (i in 1:nd) count[i] <- eval(parse(text=paste("ncid$dim$'",cdfdims[i],"'$len",sep="")))
  varsize <- count
  lon.we <- lon
  if (!is.null(x.rng)) {
    if (!silent) print(paste("Longitudes: ",min(lon[is.finite(lon)]),"-",max(lon[is.finite(lon)]),attr(lon,"unit")))
    if (!silent) print(paste("extract: ",min(x.rng),"-",max(x.rng)))
    if (min(x.rng) < 0 & max(lon > 180)) start[1] <- max(c(sum(lon < min(x.rng))),1) else  # REB fix 21.10.2005
                                         start[1] <- max(sum(lon < min(x.rng))+1)          # REB fix 21.10.2005/ REB fix 10.01.2008
### Problem indicated by Lijun Fan Jan 08 2008 - selecged area shifted one grid-box to the southwest.    
    ix <- (lon >= min(x.rng) & lon <= max(x.rng))
    if (!silent) print(paste("sum(ix)=",sum(ix)))
    lon <- lon[ix]
     attr(lon,"unit") <- eval(parse(text=paste("ncid$'dim$",cdfdims[ilon],"$units'",sep="")))
    count[1] <- length(lon)
  }
  if (!is.null(y.rng)) {
    if (!silent) print(paste("Latitudes: ",min(lat[is.finite(lat)]),"-",max(lat[is.finite(lat)]),attr(lat,"unit")))
    if (!silent) print(paste("extract: ",min(y.rng),"-",max(y.rng)))

    if (lat[1] > lat[2]) start[2] <- max(sum(lat < min(y.rng)),1) else          # REB fix 21.10.2005/ REB fix 10.01.2008
                         start[2] <- max(sum(lat < min(y.rng))+1,1)             # REB fix 21.10.2005/ REB fix 10.01.2008
    iy <- (lat >= min(y.rng) & lat <= max(y.rng))
    lat <- lat[iy]    
    attr(lat,"unit") <- eval(parse(text=paste("ncid$dim$",cdfdims[ilat],"$units",sep="")))
    count[2] <- length(lat)
  }
  if ((!is.null(z.rng)) & nd ==4  ) {
    if (lev[1] > lev[2]) start[3] <- min(sum(lev < min(z.rng)),1) else          # REB fix 21.10.2005/ REB fix 10.01.2008
                         start[3] <- min(sum(lev < min(z.rng))+1,1)             # REB fix 21.10.2005/ REB fix 10.01.2008 
    iz <- (lev >= min(z.rng) & lev <= max(z.rng))
    lev <- lev[iz]    
    count[3] <- length(lev)
    nz <- length(lev)
  }


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
        tim <- seq(tim[1],tim[1]+(nt-1)*dtim[1],by=median(dtim))
      } else {
        dt <- readline("What is the correct time step? (0 leaves tim unchanged)")
        if (dt != 0) tim <- seq(tim[1],tim[1]+nt-1,by=dt)
      }
    }
    print(paste("length(tim)=",length(tim),"nt=",nt))
  }

  if ((is.null(torg)) & (regexpr("since",t.unit)[1]>0)) {
    torg <- substr(t.unit,regexpr("since",t.unit)+6,nchar(t.unit))
    t.unit <- substr(t.unit,1,regexpr("since",t.unit)-2)
  }

  if (is.null(torg)) {
    print(paste("Time units:",t.unit," l=",min(tim[is.finite(tim)]),"-",max(tim[is.finite(tim)])))
    print("Cannot determine the time origin!")
    print("Example format: '15-Dec-1949'")
    print("NCEP reanalysis typically: 01-01-01")
    print("ERA-40 typically: 1900-01-01")
    torg <- readline("I need a time origin: ")
  } 

  if (!is.null(torg)) {
    print(paste("torg=",torg))
    yy0 <- datestr2num(torg)[1]
    mm0 <- datestr2num(torg)[2]
    dd0 <- datestr2num(torg)[3]
    if (is.na(dd0)) dd0 <- 1
    if (is.na(mm0)) mm0 <- 15  
  } else torg <- readline("Give me the time origin (format='15-Dec-1949'):")

  if (!silent) print(paste("Time origin: (year-month-day)",yy0,"-",mm0,"-",dd0))
  if ( (yy0[1]==0) & (forceBC) ) {
    if (!silent) print('There is no year zero (Press et al., Numerical recipies)')
    if (!silent) print("'> print(julday(1,1,1)-julday(1,1,-1))' gives 365")
    if (!silent) print('julday wont work unless the time is fixed')
    if (!silent) print("year0 is set to 1, and 365 days is subtracted from tim")
    if (substr(lower.case(t.unit),1,4)=="hour") {
      tim <- tim - 365*24
      t.unit <- "day"
    }
    if (substr(lower.case(t.unit),1,3)=="day") tim <- tim - 365
    if (substr(lower.case(t.unit),1,3)=="mon") tim <- tim - 12
    if (substr(lower.case(t.unit),1,4)=="year") tim <- tim - 1
    yy0 <- 1
  }

  if (!silent) print(paste("Time unit:",lower.case(t.unit)))
  if (substr(lower.case(t.unit),1,3)=="mon") {
    tim <- floor(tim)
    mm <- mod(mm0 + tim - 1,12)+1
    yy  <- yy0 + floor((tim+mm0-1)/12)
    dd <- rep(15,length(tim))
    obj.type <- "monthly.field.object"
  } else if (substr(lower.case(t.unit),1,3)=="day") {
    if (yy0!=0) mmddyy <- caldat(tim + julday(mm0,dd0,yy0)) else if (median(diff(tim)) > 29){
      year <- yy0 + floor(tim/daysayear)
      month <- mm0 + rep(1:12,ceiling(length(year)/12))[1:length(tim)] -1
      day <- rep(15,length(tim))
      mmddyy <- list(day=day,month=month,year=year)
    } else stop(paste('There is a problem with the time dimansion - I do not know what to do.',
                      'Can be fixed witrh NCO? (http://sf.net/projects/nco)'))
    mm <- mmddyy$month
    yy <- mmddyy$year
    dd <- mmddyy$day
    obj.type <- "daily.field.object"
  } else if (substr(lower.case(t.unit),1,4)=="hour") {
    if (yy0!=0) mmddyy <- caldat(tim/24 + julday(mm0,dd0,yy0)) else {
    }
    mm <- mmddyy$month
    yy <- mmddyy$year
    dd <- mmddyy$day
    tim <- tim/24
    t.unit <- "day"
    obj.type <- "daily.field.object"
  } 

  nt <- length(tim)
  if (!is.null(t.rng)) {
    if (!is.character(t.rng)) {
      start[nd] <- t.rng[1]
      if (start[nd] > nt) {
        print("Argument 't.rng':")
        print(t.rng)
        print("Detected an error: start given by 't.rng' exceeds physical record length")
        print("Are you mixing up dates (given as string argument) and index (numerical argument)?")
        stop()
      }
      tim <- tim[t.rng[1]:t.rng[2]]
      attr(tim,"unit") <- eval(parse(text=paste("ncid$dim$",cdfdims[itim],"$units",sep="")))
      count[nd] <- length(tim)
      yy <- yy[t.rng[1]:t.rng[2]]; mm <- mm[t.rng[1]:t.rng[2]]; dd <- dd[t.rng[1]:t.rng[2]]
    } else {
      yy.1<-datestr2num(t.rng[1])[1]; mm.1<-datestr2num(t.rng[1])[2]; dd.1<-datestr2num(t.rng[1])[3]
      yy.2<-datestr2num(t.rng[2])[1]; mm.2<-datestr2num(t.rng[2])[2]; dd.2<-datestr2num(t.rng[2])[3]
      it <- 1:length(tim)
      it1 <- min(it[(yy*10000 + mm*100 + dd >= yy.1*10000+mm.1*100+dd.1)],na.rm=TRUE)
      it2 <- max(it[(yy*10000 + mm*100 + dd <= yy.2*10000+mm.2*100+dd.2)],na.rm=TRUE)
      if (!silent) print(c(yy.1,mm.1,dd.1)); print(c(yy.2,mm.2,dd.2)); print(c(it1,it2))
      tim <- tim[it1:it2]
      start[nd] <- max(c(1,it1),na.rm=TRUE)
      count[nd] <- length(tim)
      yy <- yy[it1:it2]; mm <- mm[it1:it2]; dd <- dd[it1:it2]
    }       
  }
  nt <- length(tim); ny <- length(lat)
 
  if (!silent) print(cbind(start,count,varsize))
  if (!is.null(y.rng) & lat[1] > lat[length(lat)]) {
    start[2] <- varsize[2] - start[2] - count[2] + 1
  }

  start[!is.finite(start)] <- 1; start[is.element(start,0)] <- 1; start[!is.numeric(start)] <- 1
  count[!is.finite(count)] <- 1; count[is.element(count,0)] <- 1; count[!is.numeric(count)] <- 1
  if (count[1]*count[2]*count[3] > 0) {
    #if (!silent) 
    print(cbind(start,count,varsize))
    nx <- length(lon); ny <- length(lat); nt <- length(tim)
    if (!silent) print(paste("Reading",v1))
    data <- get.var.ncdf(ncid,v1,start=start,count=count)
    dim(data) <- count
#x11(); image(lon,lat,data[,,1],main="eastern H."); addland()
  } else data <- NULL
  
  if (min(x.rng) < 0 & max(lon.we > 180)) {
    if (!silent) print("read the data from western hemisphere also")
    lon.we[lon.we > 180] <- lon.we[lon.we > 180] - 360
#print(lon.we)
    start[1] <- seq(1,length(lon.we),by=1)[(lon.we >= min(x.rng)) & (lon.we < 0)][1]
    ix <- (lon.we >= min(x.rng) & lon.we < 0)
    if (sum(ix)> 0) {
      lon.we <- lon.we[ix]
#print(lon.we)
      count[1] <- length(lon.we)
      lon <- lon[!is.element(lon,lon.we)]
      if (!silent) print(cbind(start,count,varsize))
      start[!is.finite(start)] <- 1; start[is.element(start,0)] <- 1; start[!is.numeric(start)] <- 1
      count[!is.finite(count)] <- 1; count[is.element(count,0)] <- 1; count[!is.numeric(count)] <- 1
      data.w <- get.var.ncdf(ncid,v1,start=start,count=count)
# HERE 05.02.2008
      print(c(dim(data.w),NA,count))
      dim(data.w) <- count 

      print(dim( data ))
        print(dim( data.w ))
        print(start)
        print(count)

      lon <- c(lon.we,lon)
#x11(); image(lon.we,lat,data.w[,,1],main="western H."); addland()
      dat <- matrix(nrow=nt,ncol=ny*(nx+count[1]))
      dim(dat) <- c(nt,ny,nx+count[1])
      if (nd==3) {
        for (i in 1:nt) {
          if (!is.null(data)) dat[i,,] <- t(rbind(matrix(data.w[,,i],sum(ix),ny),matrix(data[,,i],nx,ny))) else
                              dat[i,,] <- t(matrix(data.w[,,i],sum(ix),ny))
        }
    } else if (nd==4) {
        for (i in 1:nt) {
          if (!is.null(data)) dat[i,,,] <- t(rbind(matrix(data.w[,,,i],sum(ix),ny),matrix(data[,,,i],nx,ny,nz))) else
                              dat[i,,,] <- t(matrix(data.w[,,,i],sum(ix),ny,nz))
        }
    }     
      attr(lon,"unit") <- eval(parse(text=paste("ncid$dim$",cdfdims[ilon],"$units",sep="")))
#print(lon); print(dim(t(dat[1,,]))); print(dim(dat)); print(length(lon)); print(length(lat))
#print(dim(data));   print(dim(data.w)); print(sum(ix))       
#x11(); image(lon,lat,t(dat[1,,]),main="east+west"); addland()
    #print("East and west combined!")
      rm(data,data.w)
    } else {
       dat <- data*NA; dim(dat) <- c(nt,ny,nx)
       for (i in 1:nt) dat[i,,] <- t(as.matrix(data[,,i]))
       rm(data)
    }
  } else {
  
 # Re-order the data: (old convention)

    if (!silent) print(c(nt,ny,nx,NA,dim(data)))
    dat <- data*NA; dim(dat) <- c(nt,ny,nx); dim(data) <- c(nx,ny,nt);  
    for (i in 1:nt) dat[i,,] <- t(as.matrix(data[,,i]))
    rm(data)
  }
  close.ncdf(ncid)
  nx <- length(lon)  
  print("dim(dat):"); print(dim(dat))

# Check for 'model dates', i.e. 360-day year
  #print("---- --- -- - Check for 'model dates', i.e. 360-day year . .. ... ...."); print(nd)
#  if (nd==3) y.test <- as.vector(dat[,1,1]) else y.test <- as.vector(dat[,1,1,1])
  y.test <- as.vector(dat[,1,1]) 
  #print(c(sum(is.finite(y.test)),length(y.test)))
  if ( (sum(is.finite(y.test)) > 100) & (daysayear != 360) ){
    ac.gcm  <- data.frame(y=y.test, x1=as.vector(cos(2*pi*tim/360)), 
                                    x2=as.vector(sin(2*pi*tim/360)))
    ac.real <- data.frame(y=y.test, x1=as.vector(cos(2*pi*tim/daysayear)), 
                                    x2=as.vector(sin(2*pi*tim/daysayear)))
    lm.gcm <- lm(y ~ x1 + x2,data=ac.gcm); r2.gcm <- summary(lm.gcm)$r.squared
    lm.real <- lm(y ~ x1 + x2,data=ac.real); r2.real <- summary(lm.real)$r.squared
    #print(summary(lm.gcm))
    #print(summary(lm.real))
    #print(force365.25) 

    if (force365.25==-1) {
      print("> > > > FORCING a '360-day' model year! < < < <")
      # BUG! REB - dating fixed 20.01.2006
#      yy <- yy0 + floor((tim +(mm0-1)*30+dd0-2)/360)        
#      mm <- mod(mm0 + floor((dd0+tim-2)/30)-1,12)+1
#      dd <- mod(dd0+tim-2,30)+1
      juldays <- caldat(tim+julday(mm0,dd0,yy0));            # REB 20.1.2006
      yy <- caldat(juldays)$year
      mm <- caldat(juldays)$month
      dd <- caldat(juldays)$day
      daysayear<- 365.25
      force365.25 <- TRUE
    }
    if (is.finite(r2.gcm) & is.finite(r2.real)) {
      if ( (r2.gcm > r2.real) & (length(rownames(table(diff(tim)))) <= 2) &
          ( (substr(lower.case(t.unit),1,3)=="day") | (substr(lower.case(t.unit),1,4)=="hour") ) & 
          !force365.25) {
        print("> > > > Detecting a '360-day' model year! < < < <")
        yy <- yy0 + floor((tim +(mm0-1)*30+dd0-2)/360)
        mm <- mod(mm0 + floor((dd0+tim-2)/30)-1,12)+1
        dd <- mod(dd0+tim-2,30)+1
        daysayear<- 360
      }
    }
  }

# Extra processing for NCEP files e.g. with Time unit: hours since 1-1-1 00:00:0.0.
  if ( ((substr(lower.case(t.unit),1,4)=="hour") |
        (substr(lower.case(t.unit),1,3)=="day")) &
        (max(diff(dd)) == 0) ) {
    if (!silent) print("Monthly data, but time unit set to 'hour'/'day'")
    if (!silent) print("Set time unit to month")
    obj.type <- "monthly.field.object"
    dd[] <- 15
    t.unit <- "month"
}

#  print("Latitude:")
  if (!is.null((attributes(lat)$"unit"))) if (attributes(lat)$"unit"=="degrees_south") lat <- lat * -1
  if (!is.null((attributes(lon)$"unit"))) if (attributes(lon)$"unit"=="degrees_west") lon <- lon * -1
  if (!is.null((attributes(lat)$"units"))) if (attributes(lat)$"units"=="degrees_south") lat <- lat * -1
  if (!is.null((attributes(lon)$"units"))) if (attributes(lon)$"units"=="degrees_west") lon <- lon * -1

  if (greenwich) {
    lon[lon > 180] <- lon[lon > 180]-360
  }

#  print("Sort longs and lats")
  x.srt <- order(lon)
  y.srt <- order(lat)
  lon <- lon[x.srt]
  lat <- lat[y.srt]
#  if (nd==3) dat <- dat[,y.srt,x.srt] else
#             dat <- dat[,,y.srt,x.srt]
  dat <- dat[,y.srt,x.srt] 
  
  nx <- length(lon)
  ny <- length(lat)
  nt <- length(tim)

  if ((max(mm) > 12) & (max(dd) <= 12)) {
    mm2 <- mm; mm <- dd; dd <- mm2; rm(mm2)
  }
  if (!silent) print(paste("First & last records:",yy[1],mm[1],dd[1],
                     "&",yy[length(yy)],mm[length(mm)],dd[length(dd)]))
  
  if (l.scale) {   
    if (!silent) print("BEFORE scale adjustment & weeding")
    if (!silent) print(summary(as.vector(dat)))
    dat[dat == miss] <- NA
    if (sum(is.na(dat))>0) print(paste(sum(is.na(dat)),"of",length(dat),
                                 " are set to 'NA'"))
      if (!silent) print("AFTER scale adjustment & weeding")
  }

  if ((l.scale) & !is.null(scal)) {
     print(paste("Scaling: dat <- dat *",scal))
     if (is.finite(scal)) dat <- dat * scal
  }
  # Have included a sanity test to detect an old 'bug': offset 273 and
  # units of deg C..
  if ( ((l.scale) & !is.null(offs))) {
      if ( (offs!=273) &
           (unit=="deg C")) {
        a <- readline(prompt="Correct an old bug? (y/n)")
        if (lower.case(a)=="y") dat <- dat + offs} else
        if (is.finite(offs)) {
          print(paste("Offset: dat <- dat +",offs))
          dat <- dat + offs
        }
  }

  if ((unit=="K") | (unit=="Kelvin") |
      (unit=="degrees Kelvin") |
      (unit=="deg K") | (unit=="degK")) {
    dat <- dat - 273
    unit <- "deg C"
  }
    if ((unit=="Pa") | (substr(strip(lower.case(unit)),1,6)=="pascal") |
      (unit=="N/m^2") |
      (unit=="N m^{-1}")) {
    dat <- dat/100
    unit <- "hPa"
  }
  if (!silent) print(summary(as.vector(dat)))

  if (!silent) print(paste("dimensions",nt,ny,nx))
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
  dat.att <- list(time.unit=t.unit, time.origin=torg, unit= unit,long.name=lon.nam,
                  filename=filename,scale.factor=scal,add.offset=offs,miss=miss,
                  daysayear=daysayear)
  attr(tim,"unit") <- t.unit
  attr(tim,"time_origin") <- torg

#print(obj.type)
#print(paste("min(diff):",min(diff(tim))))
#print(sum(is.element(diff(mm),0)))

  if ((obj.type=="daily.field.object") & (min(diff(tim))>=28) & (regular) &
      (sum(is.element(diff(mm),0))>0) & (substr(lower.case(t.unit),1,3)=="day")) {
     if (!silent) print("Problems detected in date-timeunit Quality Control")
     tim <- 1:length(tim)-1; t.unit <- "month"
     attr(tim,'unit') <- "months"
     mm <- mod(mm0 + tim - 1,12)+1
     yy <- yy0 + floor((tim+mm0-1)/12)
     dd <- rep(15,length(tim))
     obj.type <- "monthly.field.object"
     if (!silent) print("Re-setting the unit to monthly!")
     if (!silent) print(paste("New First & last dates:",yy[1],mm[1],dd[1],
                     "&",yy[length(yy)],mm[length(mm)],dd[length(dd)]))
  }

  retrieve.nc  <- list(dat=dat,lon=lon,lat=lat,tim=tim,lev=lev,
                       v.name=v.nam,id.x=id.x,id.t=id.t,
                       yy=yy,mm=mm,dd=dd,n.fld=1,
                       id.lon=rep(v.nam,nx),id.lat=rep(v.nam,ny),
                       attributes=dat.att,filename=filename)
  class(retrieve.nc) <- c("field",obj.type)
  invisible(retrieve.nc)
}



#----------------------------------------------------------------

fixField <- function(x,torg=NULL,t.unit=NULL,scal=NULL,offs=NULL, 
                     x.rng=NULL,y.rng=NULL,z.rng=NULL,t.rng=NULL,greenwich=TRUE) {
  tim <- x$tim; lon <- x$lon; lat <- x$lat
  mm <- x$mm; yy <- x$yy; dd <- x$dd
  
  if (!is.null(torg)) {
    dsh <- instring("-",torg)
    print(dsh)
    yy0 <- as.numeric(substr(torg,dsh[2]+1,dsh[2]+4))
    dd0 <- as.numeric(substr(torg,1,dsh[1]-1))
    mm0 <- switch(lower.case(substr(torg,dsh[1]+1,dsh[2]-1)),
                  "jan"=1,"feb"=2,"mar"=3,"apr"=4,"may"=5,"jun"=6,
                  "jul"=7,"aug"=8,"sep"=9,"oct"=10,"nov"=11,"dec"=12)
#    print(c(torg,dd0,mm0,yy0))

    print(paste("Time origin: (year-month-day)",yy0,"-",mm0,"-",dd0))
    if (yy0[1]==0) {
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

    if (is.null(t.unit)) t.unit <- x$dat.att$t.unit
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
  } else torg <- x$dat.att$torg
  
  if (!is.null(scal)) x$dat <- x$dat + scal

  if (!is.null(offs)) x$dat <- x$dat + offs

  if (greenwich) {
    lon[lon > 180] <- lon[lon > 180]-360
  }
  if (!is.null(x.rng)) {
    print(range(lon))
    print("Extract longitudes:")
    print(x.rng)
    x.keep <- (lon >= min(x.rng)) & (lon <= max(x.rng))
    if (nd==3) dat <- dat[,,x.keep] else
               dat <- dat[,,,x.keep] 
    lon <- lon[x.keep]
    id.lon <- id.lon[x.keep]
    id.x <- id.x[,x.keep]
    x.srt <- order(lon)
    lon <- lon[x.srt]
    if (nd==3) dat <- dat[,,x.srt] else
               dat <- dat[,,,x.srt]
  }

  if (!is.null(y.rng)) {
    print(range(lat))
    print("Extract latitudes:")
    print(y.rng)
    y.keep <- (lat >= min(y.rng)) & (lat <= max(y.rng))
    if (nd==3) dat <- dat[,y.keep,] else
               dat <- dat[,,y.keep,] 
    lat <- lat[y.keep]
    id.lat <- id.lat[y.keep]
    id.x <- id.x[y.keep,]
    y.srt <- order(lat)
    lat <- lat[y.srt]
    if (nd==3) dat <- dat[,y.srt] else
               dat <- dat[,,y.srt]
  }

  if (!is.null(t.rng)) {
    print(range(tim))
    print("Extract tims:")
    print(tim)
    t.keep <- (tim >= min(t.rng)) & (tim <= max(t.rng))
    if (nd==3) dat <- dat[t.keep,,] else
               dat <- dat[t.keep,,,] 
    tim <- tim[t.keep]
    id.t <- id.t[t.keep]
    yy <- yy[t.keep]; mm <- mm[t.keep]; dd <- dd[t.keep]
  }

  x$tim <- tim; x$lon <- lon; x$lat <- lat
  x$yy <- yy; x$mm <- mm; x$dd <- dd
  x$dat.att$t.unit <- t.unit
  x$dat.att$torg <- torg
  x$dat.att$scale.factor <- scal
  x$dat.att$add.offset <- offs
  if (is.null(x$dat.att$fixes)) x$dat.att$fixes <- paste("fixField ",date(),": ",torg,t.unit,scal,offs,sep="") else
                                x$dat.att$fixes <- paste(x$dat.att$fixes,date(),": ",torg,t.unit,scal,offs,sep="")
  invisible(x)
}
