lagStation <- function(x,lag=0) {
  nt <- length(x$val)
  dims <- dim(x$val)
  y <- as.vector(t(x$val))
  yshift <- rep(NA,nt)
  if (lag >= 0) yshift[1:(nt-lag)] <- y[(1+lag):nt] else
                yshift[(1-lag):nt] <- y[1:(nt+lag)]
  x$val <- t(matrix(yshift,dims[2],dims[1]))
  if (!is.null(attr(x,'lagStation'))) attr(x,'lagStation')<- paste(attr(x,'lagStation'),lag) else
                                      attr(x,'lagStation') <- lag
  invisible(x)
}

ds2station <- function(x,what="scenario") {
   ele <- switch(lower.case(x$v.name),
                 "t" = 101, "temp" = 101, "t2m" = 101,
                 "p" = 601, "precip" = 601, "rain" = 601, "precipitation" = 601)
   if (what=="scenario") {
     y <- station.obj(x=x$pre.gcm,yy=x$yy.gcm,mm=x$mm.gcm,obs.name=x$v.name,
                      unit=x$unit,ele=ele,lat=x$lat.loc,lon=x$lon.loc,alt=x$alt.loc,
                      location=x$location,ref="ds2station")
   } else {
     y <- station.obj(x=x$pre.cal,yy=x$yy.cal,mm=x$mm.cal,obs.name=x$v.name,
                      unit=x$unit,ele=ele,lat=x$lat.loc,lon=x$lon.loc,alt=x$alt.loc,
                      location=x$location,ref="ds2station")
 } 
}

mergeEOF <- function(eof1,eof2,plot=TRUE,silent=FALSE,method="lm",
                     match="time",cut.off=8,adjust=TRUE) {
  
  if ((class(eof1)[1]!="eof") | (class(eof2)[1]!="eof")) {
    print(class(eof1)); print(class(eof2))
    stop('Need two eof objects')
  }
   
  if (match=="space") method <- "project"
  i1 <- is.element(eof1$yy+(eof1$mm-0.5)/12+(eof1$dd-0.5)/365, 
                   eof2$yy+(eof2$mm-0.5)/12+(eof2$dd-0.5)/365)
  i2 <- is.element(eof2$yy+(eof2$mm-0.5)/12+(eof2$dd-0.5)/365, 
                   eof1$yy+(eof1$mm-0.5)/12+(eof1$dd-0.5)/365)
  if (((sum(i1)==0) | (sum(i2)==0)) & method=="lm") {
    print(c(range(eof1$yy),NA,range(eof1$mm),NA,range(eof1$dd)))
    print(c(range(eof2$yy),NA,range(eof2$mm),NA,range(eof2$dd)))
    print("Use the projection method and set argument match='space'.")
    stop('No overlapping times')
  } else {
    if (!silent) print(paste("Common points: N=",sum(i1)))
    if (!silent) print(c(range(eof2$yy[i2]),NA,range(eof2$mm[i2]),NA,range(eof2$dd[i2])))
  }
  eof.match <- eof1
  eof.match$PC <- eof1$PC[!i1,]
  X <- as.matrix(eof1$PC[i1,] %*% diag(eof1$W))

  if (method=="lm") {                                       # Method section: "lm"
    lm.str <- "lm(y ~ X.1"
    for (i in 2:length(eof1$W)) {
      lm.str <- paste(lm.str,' + X.',i,sep="")
    }
    lm.str <- paste(lm.str,',data=cal)',sep="")
    #print(lm.str)
    r.squared <- eof2$W + NA
    cut.off <- min(c(cut.off,length(eof2$W)))
    #print(paste("cut-off=",cut.off))

    for (i in 1:cut.off) {
      #print(paste("<<< i=",i,">>>"))
      y <- as.vector(eof2$PC[i2,i] * eof2$W[i])
      cal <- data.frame(y = y, X = X)
      idp <- data.frame(X = as.matrix(eof1$PC[!i1,] %*% diag(eof1$W)))
      #print(summary(cal))

      match <- eval(parse(text=lm.str))
      step.wise <- step(match,trace=0)
      stats <- summary(step.wise)

# Determine which PCs were selected in the step procedure
# Note, some of the higher modes are truncated

      c<-as.character(step.wise$call[2])
      c<-unlist(strsplit(c," \\~ "))
      c<-c[2]
      c<-paste(unlist(strsplit(c," \\+ ")),' ')
      #print(c)
      incl<- rep(FALSE,length(eof1$var.eof))
      for (iv in 1:length(eof1$var.eof)) {
        if ( (!is.na( charmatch(paste('X',as.character(iv),' ',sep=""),c) )) |
             (!is.na( charmatch(paste('X.',as.character(iv),' ',sep=""),c) )) )
        {
          incl[iv]<-TRUE
        }
      }
 }
  # print(incl)

# Note that the intercept is included in lm.coe, but not in the
# coefficients held by c.

      lm.coe <- coef(step.wise)

# Find the predictor patterns

      preds2D<-eof1$EOF
      dims <- dim(preds2D)
      if (length(dims) > 2) dim(preds2D)<-c(dims[1],dims[2]*dims[3])

      if (!silent) print(paste("Reconstruct the spatial patterns",i))
      if (!silent) print(lm.coe)
      if (length(lm.coe)<2) {
        stop('Poor match!')
      }
     if (length(step.wise$coefficients) > 1) { 
      i.last <- 0
      id <- row.names(table(eof1$id.x))
      ordr <- rep(NA,length(id))
      for (iii in 1:length(id)) {
        ordr[iii] <- min((1:length(eof1$id.x))[is.element(eof1$id.x,id[iii])])
      }
      id<-id[order(ordr)]
 
      y.hat <- predict(step.wise, newdata=idp)
      y.fit <- predict(step.wise)
      #W <- as.numeric(step.wise$coefficients[-1])*eof1$W[incl[1:length(eof1$W)]]
      W <- as.numeric(step.wise$coefficients[-1])
       for (ii in 1:eof1$n.fld) {
        print(paste("field",ii))
        i.lon <- eof1$id.lon == id[ii]
        i.lat <- eof1$id.lat == id[ii]
        ny<-eof1$size[2,ii]
        nx<-eof1$size[3,ii]
        i.fld <- seq(i.last+1,i.last+ny*nx,by=1)
        i.last <- max(i.fld)
 
        EOF.1 <- t(preds2D[,i.fld])
        EOF.1 <-  EOF.1[,incl[1:length(eof1$W)]]
 
        expr <- paste("X.",ii," <- EOF.1 %*% W + step.wise$coefficients[1]",sep="")
        eval(parse(text=expr))
 
        expr <- paste("dim(X.",ii,") <- c(ny,nx)",sep="")
        eval(parse(text=expr))
        eval(parse(text=paste("lon.",ii," <- eof1$lon[i.lon]",sep="")))
        eval(parse(text=paste("lat.",ii," <- eof1$lat[i.lat]",sep="")))
 
        r.squared[i] <- stats$r.squared
        if ((stats$r.squared < 0.45) & (cut.off > i)) {
          cut.off <- i
          print(paste(">>> i=",i,"cut.off=",cut.off," R^2=",round(r.squared[i],2)))
        }
         if (!silent) print(summary(match))
  
        eof.match$EOF[i,i.fld] <- eval(parse(text=paste("X.",ii,sep="")))
      }                                      # end ii-loop
      #print("Store predicted PC")
      eof.match$PC[,i] <- y.hat/eof2$W[i]

    } else print(paste("Bad match for field",i))
                                                                        # end method section: "lm"
  } else {

# Use projection from G. Strang (1988), "Linear algrebra and its
#     applications", Hartcourt Brace & Company, 3rd ed. (p.147):
# psi <- t(X) %/% ( t(X)%*%(X) )
# y.hat <- t(psi %*% t(X)) %*% t(Y)
    print("Use projection from G. Strang (1988):")
  
    if (match=="time") {
      X <- eof2$PC[i2,1:cut.off]
    } else {
      X <- eof2$EOF[1:cut.off,]
    }
    psi <- t(X) %/% ( t(X)%*%(X) )
    eof.match$PC <- t(psi %*% t(X)) %*% t(eof1$PC[!i1,1:cut.off]) 
    eof.match$EOF <-  t(psi %*% t(X)) %*% eof1$EOF[1:cut.off,] 
  }

# Arrange the data
print("Arrange the data")
  eof <- eof2
  print(paste("cut-off=",cut.off))
  #print(dim(eof1$EOF))
  #print(dim(eof.match$EOF))

  nt1 <- sum(!i1); nt2 <- length(eof2$id.t)
  np1<- length(eof$EOF[1,]); np2 <- length(eof.match$EOF[1,])
  eof$EOF <- cbind(matrix(eof$EOF[1:cut.off,],cut.off,np1),
                   matrix(eof.match$EOF[1:cut.off,],cut.off,np2)) 
  eof$PC <- rbind(matrix(eof.match$PC[,1:cut.off],nt1,cut.off),
                  matrix(eof2$PC[,1:cut.off],nt2,cut.off))
  eof$r.squared <- cbind(matrix(rep(r.squared[1:cut.off],nt1),cut.off,nt1),
                         matrix(rep(1,cut.off*nt2),cut.off,nt2))

  id1 <- paste(eof2$id.x[1],".1",sep="")
  id.match <- paste(eof1$id.x[1],".merge",sep="")
  eof$id.x <- c(eof2$id.x,
                paste(eof1$id.x,".merge",sep=""))
  eof$id.lon <- c(eof2$id.lon,
                  paste(eof1$id.lon,".merge",sep=""))
  eof$id.lat <- c(eof2$id.lat,
                  paste(eof1$id.lat,".merge",sep=""))
  eof$var <- eof$var[1:cut.off]
  eof$W <- eof$W[1:cut.off]
  eof$dW <- eof$dW[1:cut.off]
  eof$n.fld <- eof$n.fld + 1
  eof$lon <- c(eof2$lon,eof.match$lon)
  eof$lat <- c(eof2$lat,eof.match$lat)
  eof$size <- cbind(c(eof2$size[1:3]),c(sum(!i1),eof1$size[2:3]))
  eof$yy <- c(eof1$yy[!i1],eof2$yy)
  eof$mm <- c(eof1$mm[!i1],eof2$mm)
  eof$dd <- c(eof1$mm[!i1],eof2$dd)
  eof$tim <- c(eof1$tim[!i1],eof2$tim) 
  eof$id.t <- c(eof1$id.t[!i1],eof2$id.t)
  class(eof) <- class(eof2)
  if (adjust) eof <- adjust.eof(eof)
  
  #print(table(eof$id.x))
  #print(table(eof$id.lon))
  #print(table(eof$id.lat))
  #print(dim(eof$EOF))
  #print(eof$size)

  if (plot) {
     plotEOF(eof)
 #   newFig()
 #   plot(eof$yy + (eof$mm-0.5)/12,eof$PC[,1],lwd=2,type="l")
 #   grid()
 #   lines(eof1$yy + (eof1$mm-0.5)/12,eof1$PC[,1],lwd=2,lty=2,col="blue")
 #   lines(eof2$yy + (eof2$mm-0.5)/12,eof2$PC[,1],lwd=2,lty=2,col="red")
  }

#  print("mergeEOF: setting class of eof explicitly")
#  print(class(eof))
  invisible(eof)
}

  mapEOF <- function(x,i.eof=1,nlevs=5,add=FALSE,
                   col=c("red","blue","darkgreen","steelblue"),lwd=2,lty=1) {
  map <- map.eof(x,i.eof=i.eof,nlevs=nlevs,add=add,col=col,lwd=lwd,lty=lty)
  invisible(map)
} 


adjust.eof <- function(x) {
  rn <- row.names(table(x$id.t))
  ordr <- rep(NA,length(rn))
  for (i in 1:length(rn)) {
    ordr[i] <- min((1:length(x$id.t))[is.element(x$id.t,rn[i])])
  }
  #print(ordr)
  rn<-rn[order(ordr)]
  #print(rn)

  dims <- dim(x$PC)
  for (i in 1:dims[2]) {
    mu <- mean(x$PC[x$id.t==x$id.t[1],i],na.rm=TRUE)
    si <- sd(x$PC[x$id.t==x$id.t[1],i],na.rm=TRUE)
     for (ii in 2:length(rn)) {
      i.cal <- (x$id.t==rn[ii])
      mu.gcm <- mean(x$PC[i.cal,i],na.rm=TRUE)
      si.gcm <- sd(x$PC[i.cal,i],na.rm=TRUE)
       x$PC[x$id.t==rn[ii],i] <- (x$PC[x$id.t==rn[ii],i] -
             mu.gcm)/si.gcm * si + mu
    }
  }
  invisible(x)
}

getgiss <- function(stnr=NULL,location=NULL,lon=NULL,lat=NULL,stations=NULL,silent=FALSE) {
  if (!silent) print("Retrieving the data from URL http://www.giss.nasa.gov/")
  if (!silent) print("Please be patient")
  if (is.null(stations)) {
    if (!silent) print("Looking up station meta-data on the station on URL")
    stations<- read.fwf("http://www.giss.nasa.gov/data/update/gistemp/station_list.txt",
               width=c(9,20,18,17,5,4,2,2,4,3,5),skip=1,as.is=TRUE,comment.char = "%",header=FALSE,
               col.names=c("number","location","Country","fill1","lat","lon",
                           "type","brighness","fill2","Contry code","Brightness index"))
    stations$lon <- as.numeric(stations$lon)/10
    stations$lat <- as.numeric(stations$lat)/10
    stations$type[stations$type=="R"] <- "Rural"
    stations$type[stations$type=="S"] <- "Surburbian"
    stations$type[stations$type=="U"] <- "Urban"
    stations$type <- as.factor(stations$type)
    ivalcont <- !is.element(stations$Country,"                  ")
    valcont <- seq(1,length(stations$Country),by=1)[ivalcont]
    for (i in 1:sum(ivalcont)) {
      imatchcont <- is.element(stations$Contry.code,stations$Contry.code[valcont[i]])
      stations$Country[imatchcont] <- stations$Country[valcont[i]]
    }
  }

  if (!silent) print("Got the station meta-data!")
  if (sum(is.element(stations$number,stnr))==0) stnr <- NULL else {
    locmatch <- is.element(stations$number,stnr)
     print(paste("Found",stations$location[locmatch],"stnr=",stnr," lon=",
                 stations$lon[locmatch]," lat=",stations$lat[locmatch],
                 "  type=",as.character(stations$type[locmatch]),
                 " country=",stations$Country[locmatch]))
  }
  if (is.null(stnr) & !is.null(location)) {
     locmatch <- is.element(substr(lower.case(stations$location),2,nchar(location)+1),
                    lower.case(location))
     if (sum(locmatch)>0) {
       stnr <- stations$number[locmatch]
       print(paste("Found",stations$location[locmatch],"stnr=",stnr," lon=",
                   stations$lon[locmatch]," lat=",stations$lat[locmatch],
                   "  type=",as.character(stations$type[locmatch]),
                   " country=",stations$Country[locmatch]))
     } else {
       print(paste("Did not find '",location,"'",sep=""))
       print(substr(lower.case(stations$location),2,nchar(location)+1))
     }
     if (length(stnr)>1) {
       i <- as.numeric(readline(paste("Which of these ( 1 -",length(stnr),")? ")))
       stnr <- stnr[i]
     }
  }
  if (!is.null(lon) & !is.null(lat) & is.null(stnr)) {
     if (!silent) print(paste("Find the nearest station to ",lon,
                              "E and ",lat,"N.",sep=""))
     dist <- distAB(lon,lat,stations$lon,stations$lat)
     distmatch <- dist == min(dist,na.rm=TRUE)
     distmatch[is.na(distmatch)] <- FALSE
     if (sum(distmatch,na.rm=TRUE)>0) {
       stnr <- stations$number[distmatch]
       locmatch <- distmatch
       print(paste("Found",stations$location[locmatch],"stnr=",stnr," lon=",
                   stations$lon[locmatch]," lat=",stations$lat[locmatch],
                   "  type=",as.character(stations$type[locmatch]),
                   " country=",stations$Country[locmatch]))
     }
  }

  if (!is.null(stnr)) {
    contcode<- stations$Contry.code[is.element(stations$number,stnr)]
    lat<- stations$lat[is.element(stations$number,stnr)]
    lon<- stations$lon[is.element(stations$number,stnr)]
    location<- stations$location[is.element(stations$number,stnr)]
    country<- stations$Country[is.element(stations$number,stnr)]
    type<- stations$type[is.element(stations$number,stnr)]
    if (nchar(stnr)==8) stnr <- paste("0",stnr,sep="")
    fname<-paste("http://www.giss.nasa.gov/data/update/gistemp/TMPDIR/tmp.",
                 contcode,stnr,".1.1/",contcode,stnr,".1.1.txt",sep="")
    if (!silent) print(fname)
    data<-read.table(fname,header=TRUE)
    yy <- data[,1]
    data[abs(data) > 99] <- NA

    t2m <- station.obj(x=as.matrix(data[,2:13]),yy=yy,
                       ele=101,station=stnr,lat=lat,lon=lon,alt=NA,unit="deg C",
                       location=location,wmo.no=NA,country=country,obs.name="Temperature",
                       ref="Hansen, J. et al. 1999. J.Geophys.Res. 104, 30997-31022")
    t2m$type <- as.character(type)
    t2m$fname <- fname
    invisible(t2m)    
  } else {
    print("Did not find a station!")
    print("Returning a list of available stations instead")
    invisible(stations)
  }
  
}



