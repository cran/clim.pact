mixFields <- function(field.1,field.2,mon=NULL,
                       interval=NULL) {

  if (class(field.1)[1] != class(field.2)[1]) {
    print(class(field.1)[1])
    print(class(field.2)[1])
    stop("The objects must have the same class")
  }
  if (!is.null(interval)) {
    if (length(interval)!=2) stop('interval must be a vector of length 2')
    if ( (interval[1] > max(c(field.1$yy,field.2$yy))) |
         (interval[2] < min(c(field.1$yy,field.2$yy))) ) {
      print('Time intervals covered by the fields:')
      print(range(field.1$yy))
      print(range(field.2$yy))
      print(interval)
      stop('Check the interval: interval=c(YYYY.1, YYYY.2)')
    }
  }

  tim.unit1 <- attr(field.1$tim,"unit")
  tim.torg1 <- attr(field.1$tim,"time_origin")
  tim.unit2 <- attr(field.2$tim,"unit")
  tim.torg2 <- attr(field.2$tim,"time_origin")
  if (lower.case(substr(tim.unit1,1,3)) != lower.case(substr(tim.unit2,1,3))) {
    print("mixFields: Time units=")
    print(c(tim.unit1,tim.unit2))
    stop('The time units must match')
  }
  
  dims.1 <- dim(field.1$dat)
  if (length(dims.1)>2) dim(field.1$dat) <- c(dims.1[1],dims.1[2]*dims.1[3])
  dims.2 <- dim(field.2$dat)
  if (length(dims.2)>2) dim(field.2$dat) <- c(dims.2[1],dims.2[2]*dims.2[3])
  
  if (!is.null(mon)) {
    if (is.null(interval)) i.mm <- is.element(field.1$mm,mon) else
      i.mm <- is.element(field.1$mm,mon) & (field.1$yy >= interval[1]) &
              (field.1$yy <= interval[2])
    field.1$dat <- field.1$dat[i.mm,]
    field.1$tim <- field.1$tim[i.mm]
    field.1$yy <- field.1$yy[i.mm]
    field.1$mm <- field.1$mm[i.mm]
    field.1$dd <- field.1$dd[i.mm]
    field.1$id.t <- field.1$id.t[i.mm]
    if (is.null(interval)) i.mm <- is.element(field.2$mm,mon) else
      i.mm <- is.element(field.2$mm,mon) & (field.2$yy >= interval[1]) &
              (field.2$yy <= interval[2])    
    field.2$dat <- field.2$dat[i.mm,]
    field.2$tim <- field.2$tim[i.mm]
    field.2$yy <- field.2$yy[i.mm]
    field.2$mm <- field.2$mm[i.mm]
    field.2$dd <- field.2$dd[i.mm]
    field.2$id.t <- field.2$id.t[i.mm]
  }

  # Match the times:
  
  i1<-is.element(field.1$yy*10000+field.1$mm*100+field.1$dd,
                 field.2$yy*10000+field.2$mm*100+field.2$dd)
  i2<-is.element(field.2$yy*10000+field.2$mm*100+field.2$dd,
                 field.1$yy*10000+field.1$mm*100+field.1$dd)
  
  field.1$dat <- field.1$dat[i1,]
  field.1$tim <- field.1$tim[i1]
  field.1$yy <- field.1$yy[i1]
  field.1$mm <- field.1$mm[i1]
  field.1$dd <- field.1$dd[i1]
  field.1$id.t <- field.1$id.t[i1]
  field.2$dat <- field.2$dat[i2,]
  field.2$tim <- field.2$tim[i2]
  field.2$yy <- field.2$yy[i2]
  field.2$mm <- field.2$mm[i2]
  field.2$dd <- field.2$dd[i2]
  field.2$id.t <- field.2$id.t[i2] 
  
  nt.1 <- length(field.1$tim)
  nt.2 <- length(field.2$tim)
  nx.1 <- length(field.1$lon)
  nx.2 <- length(field.2$lon)
  ny.1 <- length(field.1$lat)
  ny.2 <- length(field.2$lat)

# This should never occur!  
  if (nt.1!=nt.2) {
    print("SOMETHING WENT VERY WRONG!")
    print(table(field.1$id.t))
    print(table(field.1$yy))
    print(table(field.1$mm))
    print(table(field.1$dd))
    print(table(field.2$id.t))
    print(table(field.2$yy))
    print(table(field.2$mm))
    print(table(field.2$dd))
    stop(paste("The fields must have the same time dimension: nt1=",
               nt.1," nt2=", nt.2))
  }

  dim(field.1$dat)<-c(nt.1,ny.1*nx.1)
  dim(field.2$dat)<-c(nt.2,ny.2*nx.2)

  field.mxf<-cbind(field.1$dat,field.2$dat)

  tim <- field.1$tim
  attr(tim,"unit") <- tim.unit1
  attr(tim,"time_origin") <- tim.torg1
  yy <- field.1$yy
  mm <- field.1$mm
  dd <- field.1$dd
  lon <- c(field.1$lon,field.2$lon)
  lat <- c(field.1$lat,field.2$lat)
  id.x <- as.vector(cbind(field.1$id.x,field.2$id.x))
  id.lon <- c(rep(field.1$id.x[1],nx.1),rep(field.2$id.x[1],nx.2))
  id.lat <- c(rep(field.1$id.x[1],ny.1),rep(field.2$id.x[1],ny.2))
  id.t <- paste(field.1$id.t,"+",field.2$id.t,sep="")
  var.name <- c(field.1$v.name,field.2$v.name)
  result  <- list(dat=field.mxf,lon=lon,lat=lat,tim=tim,v.name=var.name,
                  id.t=id.t,id.x=id.x,yy=yy,mm=mm,dd=dd,
                  n.fld=field.1$n.fld+field.2$n.fld,
                  id.lon=id.lon,id.lat=id.lat)
  class(result) <- c(class(field.1),"mix.fields")
  invisible(result)
}
