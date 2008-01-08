# R.E. Benestad, met.no, Oslo, Norway 22.05.2002
# rasmus.benestad@met.no
#-------------------------------------------------------------------
# Estimate anomalies

anomaly.station <- function(obs,period=c(1961,1990),param=c("t2m","precip")) {


cmon<-c("Jan","Feb","Mar","Apr","May","Jun",
        "Jul","Aug","Sep","Oct","Nov","Dec")

if (lower.case(class(obs)[2])=="monthly.station.record") {
  ny <- length(obs$yy)
  value <- t(obs$val)
  clim <- rep(NA,12)
  if (!is.null(period)) ii <- ((obs$yy>=min(period)) & (obs$yy<=max(period))) else
                        ii <- is.finite(obs$yy)
  for (im in 1:12) {
        clim[im] <- mean(value[im,ii],na.rm=TRUE)
        value[im,] <- value[im,] - clim[im]
      }
  obs$val <- t(value)
  obs$clim <- clim
  obs$obs.name  <-  paste(obs$obs.name,"anomaly")
  }  else if (lower.case(class(obs)[2])=="daily.station.record") {
    if (!is.null(attr(obs$tim,"daysayear"))) daysayear <- attr(obs$tim,"daysayear") else
                                             daysayear <- 365.25
    #print(paste("obs$",param[1],sep=""))
    y <- eval(parse(text=paste("obs$",param[1],sep="")))
    nt <- length(y)
    jtime <- julday(obs$mm, obs$dd, obs$yy) - julday(1,1,1950)
    #print(c(length(y),length(cos(2*pi*jtime/daysayear))))

    ac.mod <- data.frame(c1=cos(2*pi*jtime/daysayear),s1=sin(2*pi*jtime/daysayear),
                         c2=cos(4*pi*jtime/daysayear),s2=sin(4*pi*jtime/daysayear),
                         c3=cos(6*pi*jtime/daysayear),s3=sin(6*pi*jtime/daysayear))

    ac.fit <- lm(y ~ c1 + s1 + c2 + s2 + c3 + s3, data=ac.mod)
    clim <- y; clim[is.finite(y)]<-ac.fit$fit
    #print(c(dim(ac.mod),NA,length(y),length(clim),length(ac.mod$c1),length(ac.fit$fit)))
    obs$abs.1 <- y
    obs$clim.1 <- clim
    #obs$t2m <- obs$t2m - clim
    cline<-paste("obs$",param[1]," <- obs$",param[1]," - clim",sep="")
    #print(cline)
    eval(parse(text=cline))
    obs$AC.model.1 <- ac.fit
    obs$ac.mod <- ac.mod


    y <- eval(parse(text=paste("obs$",param[2],sep="")))
    if (sum(is.finite(y)) >0) {
      #print(c(length(y),length(cos(2*pi*jtime/daysayear))))
      #print(range(y,na.rm=TRUE))
      ac.fit <- glm(y ~ c1 + s1 + c2 + s2 + c3 + s3, data=ac.mod)
      clim[] <- NA; clim[is.finite(y)]<-ac.fit$fit
      obs$abs.2 <- y
      obs$clim.2 <- clim
      #obs$precip <- obs$precip - clim
      cline <- paste("obs$",param[2]," <- obs$",param[2]," - clim",sep="")
      #print(cline)
      eval(parse(text=cline))
      obs$AC.model.2 <- ac.fit
    }
  }
#print("exit")
  invisible(obs)
}

addClim.station <- function(x,param=c("t2m","precip")) {
  if ((lower.case(class(x)[2])=="monthly.station.record") & !is.null(x$clim)) {
    for (im in 1:12) {
      x$value[,im] <- x$value[,im] - x$clim[im]
    }
    x$clim[] <- 0
  } else if ((lower.case(class(x)[2])=="daily.station.record")) {
    if (!is.null(attr(x$tim,"daysayear"))) daysayear <- attr(x$tim,"daysayear") else
                                             daysayear <- 365.25
      jtime <- julday(x$mm, x$dd, x$yy) - julday(1,1,1950)
      AC.mod <- data.frame(c1=cos(2*pi*jtime/daysayear),s1=sin(2*pi*jtime/daysayear),
                           c2=cos(4*pi*jtime/daysayear),s2=sin(4*pi*jtime/daysayear),
                           c3=cos(6*pi*jtime/daysayear),s3=sin(6*pi*jtime/daysayear))
      ac.reconstr <- predict(x$AC.model.1,newdata=AC.mod)
      eval(parse(text=paste("x$",param[1]," <- x$",param[1]," + ac.reconstr",sep="")))
      ac.reconstr <- predict(x$AC.model.2,newdata=AC.mod)
      eval(parse(text=paste("x$",param[2]," <- x$",param[2]," + ac.reconstr",sep="")))
 }
 invisible(x)
}

preClim.station <- function(x,dd=NULL,mm=NULL,yy=NULL,param=1) {
  if ((lower.case(class(x)[2])=="monthly.station.record") & !is.null(x$clim)) {
    ac.reconstr <- rep(NA,length(mm))
    for (i in 1:length(dd)) {
      ac.reconstr[i] <- x$clim[mm[i]]
    }
  } else if ((lower.case(class(x)[2])=="daily.station.record")) { 
    if (is.null(x$AC.model.1) & (lower.case(class(x)[1])=="station")) x <- anomaly.station(x) else
       if (is.null(x$AC.model.1)) stop("preClim.station - Error: 'x' is not a valid station object!")
    if (!is.null(attr(x$tim,"daysayear"))) daysayear <- attr(x$tim,"daysayear") else
                                             daysayear <- 365.25
   jtime <- julday(mm, dd, yy) - julday(1,1,1950)
   if (param==1) AC.model <- x$AC.model.1 else AC.model <- x$AC.model.2
   AC.mod <- data.frame(c1=cos(2*pi*jtime/daysayear),s1=sin(2*pi*jtime/daysayear),
                        c2=cos(4*pi*jtime/daysayear),s2=sin(4*pi*jtime/daysayear),
                        c3=cos(6*pi*jtime/daysayear),s3=sin(6*pi*jtime/daysayear))
   ac.mod <- x$ac.mod
# Following line doesn't work:
# Error: variables 'c1', 's1', 'c2', 's2', 'c3', 's3' were specified differently from the fit
#   ac.reconstr <- predict(AC.model,newdata=AC.mod)

   c <- AC.model$coefficients
   ac.reconstr <- c[1] + c[2]*AC.mod$c1 + c[3]*AC.mod$s1 + c[4]*AC.mod$c2 + c[5]*AC.mod$s2 +
                         c[6]*AC.mod$c3 + c[7]*AC.mod$s3
  }    
  invisible(ac.reconstr)

}


daily2monthly.station <- function(obs,param="t2m",min.days.month=20,method="mean",na.rm=TRUE) {
  if (lower.case(class(obs)[2])!="daily.station.record") stop("daily2monthly.station: needs daily.sation.record")
  x <- eval(parse(text=paste("obs$",param,sep="")))

  cont <- names(obs)
  cont <-  cont[!is.element(cont,c("dd","mm","yy","obs.name","unit","station",
                                   "lat","lon","alt","ele","x.0E65N","y.0E65N",
                                   "found","location","wmo.no","start","yy0",
                                   "country","ref"))]
  if ( (length(cont)!=length(obs$unit)) | (length(cont)!=length(obs$ele)))
    stop("daily2monthly.station: number of elements doesn't match number of units/codes")
  unit <- obs$unit[is.element(cont,param)]
  ele <- obs$ele[is.element(cont,param)]
  years <- as.numeric(rownames(table(obs$yy))); ny <-  length(years)
  months <- 1:12
  val <- rep(NA,ny*12); dim(val) <- c(ny,12)
  for (iy in 1:ny) {
    for (im in 1:12) {
      this.month.year <-  is.element(obs$mm,im) & is.element(obs$yy,years[iy])
      if (sum(this.month.year,na.rm=TRUE) > min.days.month) {
        if (!na.rm) val[iy,im] <- eval(parse(text=paste(method,"(x[this.month.year])",sep=""))) else 
            val[iy,im] <- eval(parse(text=paste(method,"(x[this.month.year],na.rm=TRUE)",sep="")))         
      }
    }
  }
  Obs <-  station.obj(x=val,yy=years,obs.name=param,unit=unit,ele=ele, mm=NULL,
                        station=obs$station,lat=obs$lat,lon=obs$lon,alt=obs$alt,
                        location=obs$location,wmo.no=obs$wmo.no,
                        start=obs$start,yy0=obs$yy0,country=obs$station,
                        ref=paste(obs$ref,"+daily2monthly.station. method=",method,
                          "na.rm=",na.rm))
  invisible(Obs)
}

