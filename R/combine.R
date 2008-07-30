combineEOFs <- function(eof.1,eof.2) {
  if ( (class(eof.1)[1]!="eof") | (class(eof.2)[1]!="eof") ) stop('combineEOFs: input must be EOF-objects')

  EOF1 <- eof.1$EOF
  EOF2 <- eof.2$EOF
  d1 <- dim(EOF1)
  d2 <- dim(EOF2)
  #print(dim(EOF1))

  if ( (d1[2] !=d2[2]) | (eof.1$lat[1] != eof.2$lat[1]) |
       (eof.1$lat[length(eof.1$lat)] != eof.2$lat[length(eof.2$lat)]) |
       (eof.1$lon[1] != eof.2$lon[1]) |
      (eof.1$lon[length(eof.1$lon)] != eof.2$lon[length(eof.2$lon)]) ) {
    lonxy <- rep(eof.2$lon,length(eof.2$lat))
    latxy <- sort(rep(eof.2$lat,length(eof.2$lon)))
    eof2 <- matrix(rep(NA,d2[1],d1[2]),d2[1],d1[2])
    for (ii in 1:d2[1]) {
      eof2[ii,] <- interp(lonxy,latxy,EOF2[ii,],eof1$lon,eof1$lat)
    }
    EOF2 <- eof2; rm(eof2)
  }
  EOF <- rbind(EOF1,EOF2)
  #print(dim(EOF))

  i1 <- is.element(eof.1$yy*10000+eof.1$mm*100+eof.1$dd,eof.2$yy*10000+eof.2$mm*100+eof.2$dd)
  i2 <- is.element(eof.2$yy*10000+eof.2$mm*100+eof.2$dd,eof.1$yy*10000+eof.1$mm*100+eof.1$dd)
  
  W <-  c(eof.1$W,eof.2$W)
  PC <- cbind(eof.1$PC[i1,],eof.2$PC[i2,])
  dW <-  c(eof.1$dW,eof.2$dW)
  n.fld <- min(c(eof.1$n.fld,eof.2$n.fld))
  var.eof <- c(eof.1$var.eof, eof.2$var.eof)
  tot.var <- eof.1$tot.var + eof.2$tot.var
  lat <- eof.1$lat
  lon <- eof.1$lon
  l.wght <- cbind(eof.1$l.wght,eof.2$l.wght)
  id.lon <- eof.1$id.lon
  id.lat <- eof.1$id.lat
  id.x <- c(eof.1$id.x,eof.2$id.x)
  id.t <- paste(eof.1$id.t[i1],eof.2$id.t[i2])
  region <- c(eof.1$region,eof.2$region)
  yy <- eof.1$yy[i1]
  mm <- eof.1$mm[i1]
  dd <- eof.1$dd[i1]
  v.name <- c(eof.1$v.name,eof.2$v.name)
  fname <- c(eof.1$fname,eof.2$fname)
  clim <- cbind(c(eof.1$clim),c(eof.2$clim))
  size <- cbind(eof.1$size,eof.2$size)
  mon <- c(eof.1$mon,eof.2$mon)
  c.mon <- c(eof.1$c.mon,eof.2$c.mon)
  attributes <- c(eof.1$attributes,eof.2$attributes)
  srt <- order(c(1:length(eof.1$W),1:length(eof.2$W)+0.5))
  EOF <- EOF[srt,]
  PC <- PC[,srt]
  W <- W[srt]
  dW <- dW[srt]
  var.eof <- var.eof[srt]
  id <- c(eof.1$id,eof.2$id)
  id.eof <- c(rep(eof.1$v.name,length(eof.1$W)),rep(eof.2$v.name,,length(eof.2$W)))[srt]
  tim <- eof.1$tim[i1]
  
  ceof<-list(EOF=EOF,W=W,PC=PC,id=id,n.fld=n.fld,tot.var=tot.var,
          id.t=id.t,id.x=id.x,size=size,dW=dW,mon=mon,l.wght=l.wght,
          id.lon=id.lon,id.lat=id.lat,region=region,tim=tim,
          lon=lon,lat=lat,var.eof=var.eof,yy=yy,mm=mm,dd=dd,
          v.name=v.name,c.mon=c.mon,f.name=fname,clim=clim,
          attributes=attributes,id.eof=id.eof)
  class(ceof) <- c(class(eof1),"combined")
  invisible(ceof)
}

combineFields <- function(field.1,field.2) {
  if ( (class(field.1)!="field") | (class(field.2)!="field") ) stop('combineFields: input must be field-objects')
  field.mxf<-cbind(field.1$dat,field.2$dat)
  tim <- field.1$tim
  attr(tim,"unit") <- tim.unit1
  attr(tim,"time_origin") <- tim.torg1
  yy <- field.1$yy
  mm <- field.1$mm
  dd <- field.1$dd
  lon <- c(field.1$lon,field.2$lon)
  lat <- c(field.1$lat,field.2$lat)
  dim(field.1$id.x) <- c(ny.1*nx.1)
  dim(field.2$id.x) <- c(ny.2*nx.2)
  id.x <- as.vector(c(field.1$id.x,field.2$id.x))
  id.lon <- c(field.1$id.lon,field.2$id.lon)
  id.lat <- c(field.1$id.lat,field.2$id.lat)
  id.t <- paste(field.1$id.t,"+",field.2$id.t,sep="")
  var.name <- c(field.1$v.name,field.2$v.name)
  result  <- list(dat=field.mxf,lon=lon,lat=lat,tim=tim,v.name=var.name,
                  id.t=id.t,id.x=id.x,yy=yy,mm=mm,dd=dd,
                  n.fld=field.1$n.fld+field.2$n.fld,
                  id.lon=id.lon,id.lat=id.lat,attributes=field.1$attributes)
  class(result) <- c(class(field.1),"combined")  
  invisible(result)
}
