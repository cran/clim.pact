composite.field <- function(x,y,lsig.mask=TRUE,sig.lev=0.05,s=0.42,mon=NULL,
                            lty=1,col="black",lwd=1,main=NULL,sub=NULL) {
  library(ctest)

if ((class(x)[1]!="field") & (class(x)[1]!="monthly.field.object") &
    (class(x)[1]!="daily.field.object")){
  stop("x must be a 'field' object.")
}

 cmon<-c('Jan','Feb','Mar','Apr','May','Jun',
         'Jul','Aug','Sep','Oct','Nov','Dec') 
  descr <- 'Composite:'
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

if (is.null(class(y))) class(y) <- 'vector'  
if (class(y)[1]=="station") {
  y.ts <- as.vector(t(y$val))
  yy <- sort(rep(y$yy,12))
  mm <- rep(1:12,length(y$yy))
  dd <- rep(15,length(yy))
  i1<-is.element(yy*10000+mm*100+dd,
                 x$yy*10000+x$mm*100+x$dd)
  i2<-is.element(x$yy*10000+x$mm*100+x$dd,
                 yy*10000+mm*100+dd)
  i.plus <- (y.ts[i1] >= mean(y.ts[i1],na.rm=TRUE)+s*sd(y.ts[i1],na.rm=TRUE))
  i.minus<- (y.ts[i1] <= mean(y.ts[i1],na.rm=TRUE)-s*sd(y.ts[i1],na.rm=TRUE))
} else if (length(y)==length(x$tim) & !is.logical(y)) {
    i1 <- seq(1,length(x$tim),by=1)
    i2 <- seq(1,length(x$tim),by=1)
    y.ts <- rep(0,length(i2))
    i.plus <- y > 0
    i.minus <- y < 0
} else if (length(y)==length(x$tim) & is.logical(y)) {
    i1 <- seq(1,length(x$tim),by=1)
    i2 <- seq(1,length(x$tim),by=1)
    y.ts <- rep(0,length(i2))
    i.plus <- y 
    i.minus <- !y    
} else if ((min(abs(y))>= min(x$yy)) & (max(abs(y)) <= max(x$yy))) {
    i1 <- seq(1,length(x$tim),by=1)    
    i2 <- seq(1,length(x$tim),by=1)
    i.plus <- is.element(x$yy,y[y > 0])
    i.minus <- is.element(x$yy,-y[y < 0])
} else {
  stop("Sorry - don't know how to interpret y; use station or vector of years")
}
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
#  print(c(sum(i.plus),sum(i.minus)))
  for (j in 1:nj) {
    for (i in 1:ni) {
      vec1 <- x$dat[i2,j,i]
      if (sum(is.finite(vec1)) > 10) {
        yy <- x$yy[i2]
        yy.plus <- yy[i.plus]
        yy.minus <- yy[i.minus]
        plus  <- vec1[i.plus]
        minus  <- vec1[i.minus]
        map[j,i] <- mean(plus,na.rm=TRUE) - mean(minus,na.rm=TRUE)
        if ((length(plus) > 0) & (length(minus) > 0)) p.val[j,i] <- t.test(plus,minus)$p.value
      } 
    }
  }

  z.levs <- seq(-max(abs(as.vector(map)),na.rm=TRUE),
                 max(abs(as.vector(map)),na.rm=TRUE),length=41)
  my.col <- rgb(c(seq(0,1,length=20),rep(1,21)),
                c(abs(sin((0:40)*pi/40))),
                c(c(rep(1,21),seq(1,0,length=20))))
  if (lsig.mask) map[p.val > 0.05] <- NA
  if (sum(is.finite(map))==0) stop('No region with significance')
  if ( (is.null(main)) & (class(y)[1]=="station") ) {
    main <- paste(descr,attributes(x$dat)$"long_name",
                  "using",y$ele,"at",y$location)
  } else if (is.null(main)) main <- paste(descr,attributes(x$dat)$"long_name")
  if (is.null(sub)) sub <- date
  filled.contour(x$lon,x$lat,t(map),
                 col = my.col,levels=z.levs,
                 main=main,sub=sub,xlab="Longitude",ylab="Latitude")

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
  results <- list(map=t(map),lon=x$lon,lat=x$lat,tim=NULL,
                  date=date,description=descr,v.name=,x$v.name,
                  yy.plus=yy.plus,yy.mins=yy.minus)
  class(results) <- "map"
  attr(results,"long_name") <- attr(x$dat,"long_name")
  attr(results,"descr") <- "Composite map"
  invisible(results)
}
