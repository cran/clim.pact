# Extracts a variable and a subregion from a netcdf File
# Assumes a netCDF file with the structure:
# field(tim,lat,lon)
# R.E. Benestad, 23.09.2003
#
# Modified 27.04.2004 to base the IO on the ncdf package instead of
# netCDF (which is being phased out).
#

cdfextract <- function(filename,varname,x.rng=NULL,y.rng=NULL,t.rng=NULL,
                       greenwich=TRUE,x.nam="lon",y.nam="lat",t.nam="tim",
                       plot=TRUE,l.scale=TRUE) {

#  library(netCDF)
  library(ncdf)

  CDF <- cdfcont(filename)
  imatch <- (CDF$vars == varname)
  if (sum(imatch)>0) {
    lon <- NULL; lat <- NULL; tim <- NULL; lev <- NULL
    dat.att <- cdfcont(filename)
    ncid1 <- open.ncdf(filename)
    v1 <- ncid1$var[[1]]
    vars <- v1$name
    nvars <-  ncid1$nvars
    print(vars)
    dims <- names(ncid1$dim)
    vars <- names(ncid1$dim)
    n.dim <- ncid1$ndims
    d <- rep(0,nvars)
    dat <- NULL
    dat.att$unit <-v1$units

    eval(parse(text=paste("lon <- ncid1$dim$",names(ncid1$dim)[1],"$vals",sep="")))
    eval(parse(text=paste("lat <- ncid1$dim$",names(ncid1$dim)[2],"$vals",sep="")))
    eval(parse(text=paste("tim <- ncid1$dim$",names(ncid1$dim)[3],"$vals",sep="")))
    attr(lon,"unit") <- eval(parse(text=paste("ncid1$dim$",names(ncid1$dim)[1],"$units",sep="")))
    attr(lat,"unit") <- eval(parse(text=paste("ncid1$dim$",names(ncid1$dim)[2],"$units",sep="")))
    attr(tim,"time_origin") <- dat.att$torg
    if (!is.null(dat.att$time.unit)) attr(tim,"unit") <- dat.att$time.unit else 
     attr(tim,"unit") <-eval(parse(text=paste("ncid1$dim$",names(ncid1$dim)[itim],"$units",sep="")))     
    if (n.dim==4) {
      eval(parse(text=paste("lev <- ncid1$dim$",names(ncid1$dim)[3],"$vals",sep="")))
      eval(parse(text=paste("tim <- ncid1$dim$",names(ncid1$dim)[4],"$vals",sep="")))
      attr(lev,"unit") <- eval(parse(text=paste("ncid1$dim$",names(ncid1$dim)[3],"$units",sep="")))
      attr(tim,"time_origin") <- dat.att$torg
      attr(tim,"unit") <- eval(parse(text=paste("ncid1$dim$",names(ncid1$dim)[4],"$units",sep="")))
    }

    print("Time information:")
    if (!is.null(dat.att$time.origin)) {
      torg <-  dat.att$time.origin
    } else torg <- NULL
 
  t.unit <- attr(tim,"unit")
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
    if (substr(lower.case(t.unit),1,5)=="month") tim <- tim - 12
    if (substr(lower.case(t.unit),1,5)=="year") tim <- tim - 1
    yy0 <- 1
  }
  print(c(mm0,dd0,yy0))
    
  print(paste("Time unit:",lower.case(t.unit)))
  if (substr(lower.case(t.unit),1,5)=="month") {
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
    
    if (greenwich) {
      lon[lon > 180] <- lon[lon > 180]-360
    }
    isrtx <- order(lon)
    print(range(lon)); print(range(lat)); print(range(tim)); 
    
    ix <- is.finite(lon); iy <- is.finite(lat); it <- is.finite(tim);
    if (!is.null(x.rng)) {ix <- ((lon>=x.rng[1]) & (lon<=x.rng[2])); print(range(lon[ix]))}
    if (!is.null(y.rng)) {iy <- ((lat>=y.rng[1]) & (lat<=y.rng[2])); print(range(lat[iy]))}
    if (!is.null(t.rng)) {it <- ((tim>=t.rng[1]) & (tim<=t.rng[2])); print(range(tim[it]))}
    iisrtx <- order(lon[ix])
    lonx <- lon[ix][iisrtx]
    nt <- length(tim); ny <- length(lat); nx <- length(lon); nz <- length(lev)
    
    x1 <- min((1:nx)[ix]); x2 <- max((1:nx)[ix])
    y1 <- min((1:ny)[iy]); y2 <- max((1:ny)[iy])
    t1 <- min((1:nt)[it]); t2 <- max((1:nt)[it])
    print(paste("netCDF dimensions:",nx,ny,nt))
    print(paste("Indeces start: ",x1,y1,t1,"& stop:",x2,y2,t2))
    print(paste("Coordinates start: ",lon[x1],lat[y1],tim[t1],
                "& stop:",lon[x2],lat[y2],tim[t2]))
    if (nz==0) {
      print(paste("Data size:",sum(it),sum(iy),sum(ix)))
      dat <- rep(0,(sum(it))*sum(iy)*sum(ix))
      dim(dat) <- c(sum(it),sum(iy),sum(ix))
      print(paste("Please be patient while reading data map by map..."))
      iit <- 0
      for (it in t1:t2) {
        iit <- iit+1
        datIN <- get.var.ncdf(ncid1,v1,,start=c(1,1,it),count=c(nx,ny,1))
        #print(c(dim(datIN),NA,nx,ny))
        dim(datIN) <- c(nx,ny)
        dat[iit,,] <- t(as.matrix(datIN[ix,iy]))
        if (plot) {
          image(lon[isrtx],lat,datIN[isrtx,],
                main=paste(yy[it],"-",mm[it],"-",dd[it]),sub=filename)
          contour(lonx,lat[iy],t(dat[iit,,iisrtx]),add=TRUE)
          lines(c(min(lonx),rep(max(lonx),2),rep(min(lonx),2)),
                c(rep(max(lat[iy]),2),rep(min(lat[iy]),2),max(lat[iy])),
                lwd=3,lty=2,col="grey95")
          addland()
          grid()
        }
       }
    } else {
       print("cdfextract does not yet know how to handle 4 dimensions ... Sorry!")
       return()
    }   
    close.ncdf(ncid1)

#print("<-------- Scale and determine time stamp...")
#print(dat.att)

    if ((l.scale) & !is.null(dat.att$scale.factor)) {
    dat <- dat * dat.att$scale.factor
  }
  # Have included a sanity test to detect an old 'bug': offset 273 and
  # units of deg C..
print("Old bug fix")

  if ( ((l.scale) & !is.null(dat.att$add.offset))) {
      if ( (dat.att$add.offset!=273) &
           (dat.att$unit=="deg C")) {
        a <- readline(prompt="Correct an old bug? (y/n)")
        if (lower.case(a)=="y") dat <- dat + dat.att$add.offset} else
        dat <- dat + dat.att$add.offset
  }

print("Set ID")

  eos <- nchar(varname)
  if (instring("-",varname)> 0) {
    eos <- instring("-",varname)-1
  } else if (instring("_",varname)> 0) {
    eos <- instring("_",varname)-1
  }
  varname <- substr(varname,1,eos)
  slash <- instring("/",filename)
  dot <- instring(".",filename)
  lon <- lon[ix]; lat <- lat[iy]; tim <- tim[t1:t2]
  isrtx <- order(lon)
  lon <- lon[isrtx]; dat <- dat[,,isrtx]
  nt <- length(tim); ny <- length(lat); nx <- length(lon); nz <- length(lev)
  id.x <- matrix(rep(varname,ny*nx),ny,nx)
  id.t <- rep(substr(filename,slash[length(slash)]+1,
                     dot[length(dot)]-1),nt)              
    
###    
print("Set coordinates")

    if (length(tim[t1:t2])>1) {
      results  <- list(dat=dat,lon=lon,lat=lat,tim=tim,lev=lev,
                       v.name=varname,id.x=id.x,id.t=id.t,
                       yy=yy[t1:t2],mm=mm[t1:t2],dd=dd[t1:t2],n.fld=1,
                       id.lon=rep(varname,nx),id.lat=rep(varname,ny),
                       attributes=dat.att)
      class(results) <- c("field",obj.type)
    } else if (length(tim[t1:t2])==1) {
      results  <- list(map=t(dat),lon=lon,lat=lat,tim=tim,
                       date=paste(yy[it],"-",mm[it],"-",dd[it]),
                       description=filename,v.name=varname)
      class(results) <- "map"
    }
#    print("Saving the extracted data in cdfextract.nc")
#    r2cdf("cdfextract.nc",results)
    invisible(results)
  } else {
    if (!file.exists(filename)) print("Cannot find filename") else {
      print(paste("Cannot find",varname,"in",filename))
      print("The netCDF file contains the following:")
      print(CDF$vars)
    }
  }
}
