monthly <- function(x,param="t2m",method="mean") {
  if (class(x)[2]!="daily.field.object" & class(x)[2]!="daily.station.record") {
      print("class(x) gives:")
      print(class(x))
      stop("Need a 'daily.field.object' or a 'daily.station.record'")
    }
  mm1 <- x$mm[1]
  mm2 <- x$mm[length(x$mm)]
  years <- as.numeric(rownames(table(x$yy)))
  nt <- length(years)*12 - mm1 + mm2 -11
  yy <- sort(rep(years,12)); yy <- yy[mm1:(length(yy) + mm2 - 12)]
  mm <- mod(mm1 + 1:nt - 2,12)+1
  dd <- rep(15,nt)
  #print(rbind(yy,mm,dd))
  #print(c(range(yy),NA,range(mm),NA,range(dd),NA,nt))
  ndays <- rep(0,nt)

  if (class(x)[2]=="daily.field.object") {
# Lines particular for a daily.field.object:
    ny <- length(x$lat); nx <- length(x$lon)
    #print(c(dim(x$dat),NA,length(x$tim),ny,nx))
    dim(x$dat) <- c(length(x$tim),ny*nx)
    dat <- x$dat[1:nt,]*0
    for (it in 1:nt) {
      ii <- is.element(x$yy,yy[it]) & is.element(x$mm,mm[it])
      ndays[it] <- sum(ii,is.na=TRUE)
      #print(c(dim(x$dat[ii,]), NA,ndays[it],NA,yy[it],mm[it],NA,it))
#      dat[it,] <- colMeans(x$dat[ii,],na.rm=TRUE)
      dat[it,] <-  apply(x$dat[ii, ], MARGIN=2,FUN=method)
    }
    dim(dat) <- c(nt,ny,nx)
    x$dat <- dat[ndays > 25,,]
    x$yy <- yy[ndays > 25]; x$mm <- mm[ndays > 25]; x$dd <- dd[ndays > 25]
    nt <- length(x$yy); x$tim <- 1:nt
    x$ndays <- ndays
    x$id.t <- rep(x$id.t[1],nt)

    if (is.null(x$dat.att$fixes)) x$dat.att$fixes <- paste("monthly",date()) else
                                x$dat.att$fixes <- paste(x$dat.att$fixes,"monthly",date())
    obj.type <- "monthly.field.object"
    attr(x$tim,"unit") <- "month"
    class(x) <- c("field",obj.type)

  } else {
# Lines particular for a daily.station.object:
   dat <- rep(NA,nt); ndays <- rep(0,nt)
   for (it in 1:nt){
     ii <- is.element(x$mm,mm[it]) & is.element(x$yy,yy[it])
     ndays[it] <- sum(ii,na.rm=TRUE)
     if (ndays[it] > 25) {
        dat[it] <- round(mean(as.numeric(eval(parse(text=paste("x$",param,"[ii]",sep="")))),na.rm=TRUE),1)
     }
    #print(c(it,mm[it],yy[it],ndays[it],dat[it]))
   }
   #print(c(length(dat),length(yy),length(mm),length(ndays)))
   #print(table(ndays))
   dat <- dat[ndays > 25]; yy <- yy[ndays > 25]; mm <- mm[ndays > 25]
   #print(c(length(dat),length(yy),length(mm),length(ndays)))
   x <- station.obj(x=dat,yy=yy,mm=mm,obs.name=param,unit=x$unit,
                    station= x$station,lat=x$lat,lon=x$lon,alt=x$alt,
                    location=x$location,country=x$country,ele=x$ele)
   x$ndays <- ndays
  }
  invisible(x)
}
