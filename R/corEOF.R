corEOF <- function(x,y,lsig.mask=TRUE,sig.lev=0.05,neofs=20,
                      lty=1,col="black",lwd=1) {
#  library(ctest)

if (class(x)[1]!="eof") {
  stop("x must be an 'eof' object.")
}

if (class(y)[1]!="station") {
  stop(paste("y must be a 'monthly.station.record'",
             "object - Use  station.obj()"))
}

  descr <- 'Correlation:'
  date <- "(EOF)"
  y.ts <- as.vector(t(y$val))
  yy <- sort(rep(y$yy,12))
  mm <- rep(1:12,length(y$yy))
  dd <- rep(15,length(yy))
  i1<-is.element(yy*10000+mm*100+dd,
                 x$yy*10000+x$mm*100+x$dd)
  i2<-is.element(x$yy*10000+x$mm*100+x$dd,
                 yy*10000+mm*100+dd)
#  print(range(yy))
#  print(range(x$yy))
#  print(range(mm))
#  print(range(x$mm))
#  print(range(dd))
#  print(range(x$dd))
#  print(c(sum(i1),sum(i2)))
  ni <- length(x$lon)
  nj <- length(x$lat)
  U <- x$EOF
  dims <- dim(U)
  if (length(dims)==3) dim(U) <- c(dims[1],dims[2]*dims[3])
  neofs <- min(c(neofs,dims[1]))
#  print('UU <- t(U %*% U)')
#  print(dim(U))
#  UU <- U^2
#  WW <- x$W^2

  ya <- (y.ts[i1] - mean(y.ts[i1],na.rm=TRUE))/sd(y.ts[i1],na.rm=TRUE)

#  Va <- x$PC[i2,]
#  for (i in 1:neofs) Va[,i] <- (Va[,i] - mean(Va[,i],na.rm=TRUE))/
#                                 sd(Va[,i],na.rm=TRUE)
#  print('WWUU')
#   if (length(dims)==2) WWUU<-matrix(rep(0,neofs*dims[2]),dims[2],neofs) else
#             WWUU<-matrix(rep(0,dims[2]*dims[3]),dims[2]*dims[3],neofs)
#  WWUU.tsum <- rep(NA,neofs)
#  for (i in 1:neofs) {
#    print(c(length(WWUU[i,]),length(UU[,i]*WW[i])))
#    
#    WWUU[i,] <- UU[,i]*WW[i]
#    WWUU.tsum[i] <- sum(WWUU[i,],na.rm=TRUE)
#  }
#  denom <- sqrt(WWUU.tsum * sum(ya^2,na.rm=TRUE))
#  denom <- sqrt(sum(x$W^2) * sum(ya^2,na.rm=TRUE))
#  denom <- sqrt(x$tot.var * sum(ya^2,na.rm=TRUE))

 # Calculate the map:
  rmap <- matrix(rep(0,ni*nj),nj,ni)
  p.val <- matrix(rep(NA,ni*nj),nj,ni)
#  for (i in 1:neofs) {
#    map <- U[i,]*x$W[i]*sum(Va[,i]*ya,na.rm=TRUE)/denom[i]
#    dim(map) <- c(nj,ni)
#    rmap[is.finite(map)] <- rmap[is.finite(map)] + map[is.finite(map)]
#  }
  X.re <- t(t(U) %*% diag(x$W) %*% t(x$PC[i2,]))
  nt <- sum(i2)
#  print(sum(is.finite(ya)))
  dim(X.re) <- c(nt,nj,ni)
    for (j in 1:nj) {
      for (i in 1:ni) {
        if (sum(is.finite(X.re[,j,i]))>25) {
          a <- cor.test(X.re[,j,i],ya)
          rmap[j,i] <- as.numeric(a$estimate)
          p.val[j,i] <- as.numeric(a$p.value)
        } else {
          rmap[j,i] <- NA
          p.val[j,i] <-NA
        }
      }
    }

  # Plot
  z.levs <- seq(-max(abs(as.vector(rmap)),na.rm=TRUE),
                 max(abs(as.vector(rmap)),na.rm=TRUE),length=41)
#  print(range(z.levs))
#  print(range(x$lon))
#  print(range(x$lat))
#  print(denom)
  my.col <- rgb(c(seq(0,1,length=20),rep(1,21)),
                c(abs(sin((0:40)*pi/40))),
                c(c(rep(1,21),seq(1,0,length=20))))
  if (lsig.mask) rmap[p.val > 0.05] <- NA


  filled.contour(x$lon,x$lat,t(rmap),
                 col = my.col,levels=z.levs,
                 main=paste(descr,attributes(x$dat)$"long_name","&",
                            y$ele,"at",y$location),
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
  contour(x$lon,x$lat,t(rmap),add=TRUE,col=col,lwd=lwd,lty=lty)
  addland()
  points(y$lon,y$lat,pch=20,col="white",cex=1.2)
  points(y$lon,y$lat,pch=20,col="black",cex=0.9) 
  results <- list(map=t(rmap),lon=x$lon,lat=x$lat,tim=x$tim,
                  date=date,description=descr)
  class(results) <- "map"
  attr(results,"long_name") <- attr(x$dat,"long_name")
  attr(results,"descr") <- "Correlation map"
  invisible(results)
}
