corField <- function(x,y,lsig.mask=TRUE,sig.lev=0.05,mon=NULL,
                     lty=1,col="black",lwd=1,main=NULL,z.levs=NULL,my.col=NULL) {
  library(ctest)
  library(akima)
  
if ((class(x)[1]!="field") & (class(x)[1]!="monthly.field.object") &
    (class(x)[1]!="daily.field.object")){
  stop("x must be a 'field' object.")
}

if ((class(y)[1]!="station") & (class(y)[1]!="field") &
    (class(y)[1]!="monthly.field.object")) {
  stop(paste("y must be a 'monthly.station.record' or a 'field'",
             "object - Use  station.obj()"))
}

 cmon<-c('Jan','Feb','Mar','Apr','May','Jun',
         'Jul','Aug','Sep','Oct','Nov','Dec') 
  descr <- 'Correlation:'
  date <- ""
  if (!is.null(mon)) {
    im <- x$mm== mon
    x$dat <- x$dat[im,,]
    x$yy <- x$yy[im]
    x$mm <- x$mm[im]
    x$dd <- x$dd[im]
    x$tim <- x$tim[im]
    x$id.t <- x$id.t[im]
    date <- cmon[mon]
  }

  if (class(y)[1]=="station") y.ts <- as.vector(t(y$val)) else
  {
    l.diffgrid <- TRUE
    if ( (length(x$lon)==length(y$lon)) & (length(x$lat)==length(y$lat)) ) {
      if ( (sum(x$lon==y$lon)==length(x$lon)) &
           (sum(x$lat==y$lat)==length(x$lat)) ) l.diffgrid <- FALSE
    }
    y.ts <- rep(NA,length(y$tim)*length(x$lat)*length(x$lon))
    dim(y.ts) <- c(length(y$tim),length(x$lat),length(x$lon))
    if (l.diffgrid) {
      print("Interpolate 2nd field - please be patient :-)")
      lat.x<-rep(x$lat,length(x$lon))
      lon.x<-sort(rep(x$lon,length(x$lat)))
      for (i in 1:length(y$tim)) {
         map <- as.matrix(y$dat[i,,])
         y.ts[i,,] <- interp(lat.x,lon.x,map,y$lat,y$lon)$z
      }
    } else {
      for (i in 1:length(y$tim)) {
         y.ts[i,,] <- as.matrix(y$dat[i,,])
      }
    }
  }
    
  if (class(y)[1]=="station") {
    yy <- sort(rep(y$yy,12))
    mm <- rep(1:12,length(y$yy))
    dd <- rep(15,length(yy))
  } else {
    yy <- y$yy; mm <- y$mm;  dd <- y$dd
  }
  
  i1<-is.element(yy*10000+mm*100+dd,
                 x$yy*10000+x$mm*100+x$dd)
  i2<-is.element(x$yy*10000+x$mm*100+x$dd,
                 yy*10000+mm*100+dd)
  ni <- length(x$lon)
  nj <- length(x$lat)
  map <- matrix(rep(NA,ni*nj),nj,ni)
  p.val <- matrix(rep(NA,ni*nj),nj,ni)
#  print(range(yy))
#  print(range(x$yy))
#  print(range(mm))
#  print(range(x$mm))
#  print(range(dd))
#  print(range(x$dd))
#  print(c(sum(i1),sum(i2)))
#  print(dim(x$dat[i2,,]))
#  print(dim(y.ts[i1,,]))
  for (j in 1:nj) {
    for (i in 1:ni) {
      if (class(y)[1]=="station") r.test <- cor.test(x$dat[i2,j,i],y.ts[i1]) else
                                  r.test <- cor.test(x$dat[i2,j,i],y.ts[i1,j,i])
      map[j,i] <- r.test$estimate
      p.val[j,i] <- r.test$p.value
    }
  }

  if (is.null(z.levs)) {
    z.levs <- seq(-max(abs(as.vector(map)),na.rm=TRUE),
                   max(abs(as.vector(map)),na.rm=TRUE),length=41)
  }
  if (is.null(my.col)) {
    my.col <- rgb(c(seq(0,1,length=20),rep(1,21)),
                  c(abs(sin((0:40)*pi/40))),
                  c(c(rep(1,21),seq(1,0,length=20))))
  }
  if (lsig.mask) map[p.val > 0.05] <- NA
  
#  print(dim(map))
#  print(c(length(x$lon),length(x$lat)))

  if (is.null(attributes(x$dat)$"long_name")) attr(x$dat,"long_name") <- x$v.name
  if (is.null(attributes(y$dat)$"long_name")) attr(y$dat,"long_name") <- y$v.name
  
  if (is.null(main)) {
    if (class(y)[1]=="station") main <- paste(descr,attributes(x$dat)$"long_name","&",
                                              y$ele,"at",y$location) else
                                main <- paste(descr,attributes(x$dat)$"long_name","&",
                                              attributes(y$dat)$"long_name")
  }
  filled.contour(x$lon,x$lat,t(map),
                 col = my.col,levels=z.levs,
                 main=main,
                 sub=date,xlab="Longitude",ylab="Latitude")

# From filled.contour in base
  mar.orig <- (par.orig <- par(c("mar","las","mfrow")))$mar
  on.exit(par(par.orig))

  w <- (3 + mar.orig[2]) * par('csi') * 2.54
  layout(matrix(c(2, 1), nc=2), widths=c(1, lcm(w)))
    
  par(las = 1)
  mar <- mar.orig
  mar[4] <- 1
  par(mar=mar)
  contour(x$lon,x$lat,t(map),add=TRUE,col=col,lwd=lwd,lty=lty)
  addland()
  if (class(y)[1]=="station") {
    points(y$lon,y$lat,pch=20,col="white",cex=1.2)
    points(y$lon,y$lat,pch=20,col="black",cex=0.9)
  }
  results <- list(map=t(map),lon=x$lon,lat=x$lat,tim=x$tim,
                  date=date,description=descr)
  class(results) <- "map"
  attr(results,"long_name") <- attr(x$dat,"long_name")
  attr(results,"descr") <- "Correlation map"
  invisible(results)
}
