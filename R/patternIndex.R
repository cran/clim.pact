# R.E. Benestad, met.no, Oslo, Norway 27.08.2003
# rasmus.benestad@met.no
#------------------------------------------------------------------------

patternIndex <- function(map,field,anomaly=TRUE) {

  library(akima)
  
  size.map <- dim(map$map)
  size.field <- dim(field$dat)
  if ( (size.map[2] != size.field[2]) |
       (size.map[1] != size.field[3]) ) {
    print("Interpolate: different grids!")
    lat.x<-rep(field$lat,length(field$lon))
    lon.x<-sort(rep(field$lon,length(field$lat)))
    Z.in<-t(as.matrix(map$map))
    Z.out<-interp(lat.x,lon.x,Z.in,map$lat,map$lon)
  } else Z.out <- t(as.matrix(map$map))
  nx <- size.field[3]; ny <- size.field[2]; nt <- size.field[1]
  
  ind <- rep(NA,nt)
  X <- as.vector(Z.out)
  Y <- field$dat
  dim(Y) <- c(nt,ny*nx)
  
  if (anomaly) {
    if ( (field$mm[2]-field$mm[1]>=1) & (field$dd[2]==field$dd[1]) ) {
      print("Months")
      for (l in as.numeric(rownames(table(field$mm)))) {
        it <- is.element(field$mm,l)
        clim <- rep(colMeans(Y[it,],na.rm=TRUE),sum(it))
#        dim(clim) <- c(sum(it),ny*nx)
        dim(clim) <- c(ny*nx,sum(it))
        clim <- t(clim)
#        plot(clim[1,],type="l",lwd=2)
#        print(dim(clim))
#        print(dim(Y[it,]))
#        points(Y[it,][1,],pch=20,col="red")
        Y[it,] <- Y[it,] - clim
#        points(Y[it,][1,],pch=21,col="blue")    # Testing OK...
      }
    } else {
      print("Other time units")
      ac.mod<-matrix(rep(NA,nt*6),nt,6)
      if (substr(lower.case(attributes(field$tim)$units),1,3)=="day") jtime <- field$tim
      if (substr(lower.case(attributes(field$tim)$units),1,4)=="hour") jtime <- field$tim/24
      ac.mod[,1]<-cos(2*pi*jtime/365.25); ac.mod[,2]<-sin(2*pi*jtime/365.25)
      ac.mod[,3]<-cos(4*pi*jtime/365.25); ac.mod[,4]<-sin(4*pi*jtime/365.25)
      ac.mod[,5]<-cos(6*pi*jtime/365.25); ac.mod[,6]<-sin(6*pi*jtime/365.25)               
      for (ip in seq(1,ny*nx,by=1)) {
        if (sum(is.finite(Y[,ip])) > 0) {
          ac.fit<-lm(Y[,ip] ~ ac.mod)
          Y[,ip]<-ac.fit$residual
        } else Y[,ip]<- NA
      }
    }
  }

  for (i in 1:nt) {
    good <- is.finite(X) & is.finite(Y[i,])
    ind[i] <- cor(as.vector(Y[i,good]),X[good])
  }

  newFig()
  dx <- (max(field$yy)-min(field$yy))/200
  plot(range(field$yy + field$mm/12 + field$dd/365.25),
       c(-1,1),type="n",lwd=3,main="Pattern Index",
       sub=field$v.name,xlab="Time",ylab="Pattern Index")
  for (ii in 1:sum(iext)) text(field$yy[iext][ii],0.975*ind[iext][ii]/abs(ind[iext][ii]),
                               as.character(field$yy[iext][ii]),cex=0.6,col="grey20",srt=90)

  lines(field$yy + field$mm/12 + field$dd/365.25+dx,ind+0.01,col="grey80",lwd=2)
  lines(field$yy + field$mm/12 + field$dd/365.25,ind,lwd=2)
  lines(c(min(field$yy),max(field$yy)+1),rep(2*sd(ind,na.rm=TRUE),2),
        lty=2,col="grey50",lwd=2)
  lines(c(min(field$yy),max(field$yy)+1),rep(-2*sd(ind,na.rm=TRUE),2),
        lty=2,col="grey50",lwd=2)
  grid()
  iext <- abs(ind) > 2* sd(ind,na.rm=TRUE)
  
  pInd <- list(index=ind,yy=field$yy,mm=field$mm,dd=field$dd,map=map)
  invisible(pInd)
}
