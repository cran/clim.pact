# R.E. Benestad, met.no, Oslo, Norway 22.05.2002
# rasmus.benestad@met.no
#-------------------------------------------------------------------
# PLot data from NORDKLIMstations.

plotStation <- function(obs,l.anom=TRUE,mon=NULL,
                        leps=FALSE,out.dir="output",what="b",
                        type="l",pch=26,col="black",lwd=3,lty=3,add=FALSE) {

if (class(obs)[2]!="monthly.station.record") {
  stop(paste("The predictand must be a 'monthly.station.record'",
             "object - Use  station.obj()"))
}

if ((!obs$found) | (sum(is.finite(obs$val))==0)) stop("No valid data!")

cmon<-c("Jan","Feb","Mar","Apr","May","Jun",
        "Jul","Aug","Sep","Oct","Nov","Dec")

if (!is.null(mon)) {
  if (((length(mon)== 1)) & (mon>0)) season <- cmon[mon]
  if (((length(mon)== 1)) & (mon==0)) season <- ""
  if ((length(mon)> 1)) season <- paste(cmon[mon[1]],'-',
                                      cmon[mon[length(mon)]],sep="")
} else {
  season <-"Dec-Jan"
}

  loc <- obs$location

  
  if (is.null(mon)) {
    ny <- length(obs$yy)
    mm <- rep(1:12,ny)
    value <- t(obs$val)
    if (l.anom) {
      for (im in 1:12) {
        value[im,] <- value[im,] - mean(value[im,],na.rm=TRUE)
        
        if ((!is.null(obs$alt)) & (!is.null(obs$lon)) & (!is.null(obs$lat))) {
          sub.tit <- paste("Anomaly",round(obs$alt,2),"m a.sl.",
                         round(obs$lon,2),"degE",
                         round(obs$lat,2),"degN")
        } else sub.tit <- paste("Anomaly:",loc)
          
      }
    } else {
      if ((!is.null(obs$alt)) & (!is.null(obs$lon)) & (!is.null(obs$lat))) {
         sub.tit <- paste("Absolute",round(obs$alt,2),"m a.sl.",
                          round(obs$lon,2),"degE",
                          round(obs$lat,2),"degN")
       } else sub.tit <- paste("Absolute:",loc)
    }

    dims <- dim(value)
    dim(value) <- c(dims[1]*dims[2],1)
    yy <- sort(rep(obs$yy,12)) + (rep(seq(1,12,by=1),ny)-0.5)/12
  } else {
    yy <- obs$yy
    ny <- length(obs$yy)
    mm <- rep(mon[1],ny)
    value <- obs$val[,mon[1]]
    if (mon==c(12,1,2)) {
      obs$val[2:ny,12] <- obs$val[1:(ny-1),12]
      obs$val[1,12] <- NA
    }
    for (i in 1:ny) value[i] <- mean(obs$val[i,mon])
    
    if (is.element(obs$ele,c(101,111,121,401,601,701,801,911)))
          for (i in 1:ny) value[i] <- mean(obs$val[i,mon],na.rm=TRUE)
    if (is.element(obs$ele,c(112,602)))
          for (i in 1:ny) value[i] <- max(obs$val[i,mon],na.rm=TRUE)
    if (is.element(obs$ele,c(122)))
          for (i in 1:ny) value[i] <- max(obs$val[i,mon],na.rm=TRUE)
    
    if ((!is.null(obs$alt)) & (!is.null(obs$lon)) & (!is.null(obs$lat))) {
       sub.tit <- paste(season," - ",round(obs$alt,2),"m a.sl.",
                     round(obs$lon,2),"degE",
                     round(obs$lat,2),"degN")
    } else sub.tit <- paste(season,": ",loc,sep="")
  }

  # Polinomial trend

  y <- value
  x <- yy
  lm.tr.p<-lm(y ~ x + I(x^2) +I(x^3) + I(x^4) + I(x^5))
  pre.p.fit<-predict(lm.tr.p,data=x)
  coef.p.fit<-lm.tr.p$coefficients
  coef.p.fit[is.na(coef.p.fit)] <- 0
  der.p.fit<-c(coef.p.fit[2],2*coef.p.fit[3],3*coef.p.fit[4],
             4*coef.p.fit[5],5*coef.p.fit[6])
  tr.est.p.fit<-(der.p.fit[1] + der.p.fit[2]*yy + der.p.fit[3]*yy^2 +
                 der.p.fit[4]*yy^3 + der.p.fit[5]*yy^4)*10


  good <- !is.na(value)
  yy <- yy[good]
  pre.p.fit <- pre.p.fit[good]
  value <- value[good]

  if (!leps) {
    
#  par(ask=TRUE)
    if ((what=="t") | (what=="b")) {
      newFig()
      par(cex.sub=0.8)
      if (!add) plot(yy,value,type="l",lwd=lwd,col=col,pch=pch,lty=lty,
                     main=paste(obs$location,obs$obs.name),
                     sub=sub.tit,xlab="Time",ylab=obs$unit) else
                lines(yy,value,lwd=lwd,col=col,pch=pch,lty=lty)
      lines(yy,pre.p.fit,col="red") 
      lines(c(min(yy),max(yy)),rep(mean(value,na.rm=TRUE)+
                                   1.96*sd(value,na.rm=TRUE),2),
                                   lty=2,col="grey")
      lines(c(min(yy),max(yy)),rep(mean(value,na.rm=TRUE)-
                                   1.96*sd(value,na.rm=TRUE),2),
                                   lty=2,col="grey")
      grid()
    }

    if ((what=="d") | (what=="b")) {
      newFig()
      par(cex.sub=0.8)
      histo <- hist(value[!is.na(value)],breaks=15,lwd=3,freq=FALSE,
         main=paste(obs$location,obs$obs.name),
         sub=paste(min(round(yy,2)),"--",max(round(yy,2)),
           ":",sub.tit,xlab=obs$unit))

      x.dist <- seq(min(histo$mids),max(histo$mids),length=101)
      y.dist <- dnorm(x.dist,
                      mean=mean(value,na.rm=TRUE),
                      sd=sd(value,na.rm=TRUE))
      lines(x.dist,y.dist,col="red")
      lines(x.dist,dgamma(x.dist-min(x.dist),
            shape=mean((value-min(x.dist))^2,na.rm=TRUE)/sd(value^2,na.rm=TRUE),
            scale=sd(value^2,na.rm=TRUE)/mean(value-min(x.dist),na.rm=TRUE)),
            col="blue",lty=3)
      grid()
    } else  {
      histo <- hist(value[!is.na(value)],breaks=15,lwd=3,freq=FALSE,plot=FALSE)
      x.dist <- seq(min(histo$mids),max(histo$mids),length=101)
      y.dist <- dnorm(x.dist,
                      mean=mean(value,na.rm=TRUE),
                      sd=sd(value,na.rm=TRUE))
    }

  } else  {
    
    figname1 <- paste(obs$location,'_',abbreviate(obs$obs.name),
                      '_',season,'1.eps',sep="")
    figname2 <- paste(obs$location,'_',abbreviate(obs$obs.name),
                      '_',season,'2.eps',sep="")
    postscript(file = figname1,onefile=TRUE,horizontal=FALSE,paper="a4")
    par(ps=14,cex.sub=0.8)
    plot(yy,value,type="l",lwd=3,
         main=paste(obs$location,obs$obs.name),
         sub=sub.tit,xlab="Time",ylab=obs$unit)
    lines(yy[!is.na(y)],pre.p.fit,col="red")
    lines(c(min(yy),max(yy)),rep(mean(value,na.rm=TRUE)+
                                 1.96*sd(value,na.rm=TRUE),2),
          lty=2,col="grey")
    lines(c(min(yy),max(yy)),rep(mean(value,na.rm=TRUE)-
                                 1.96*sd(value,na.rm=TRUE),2),
          lty=2,col="grey")
    grid()
    dev.off()

    postscript(file = figname2,onefile=TRUE,horizontal=FALSE,paper="a4")
    par(ps=14,cex.sub=0.8)
    histo <- hist(value,breaks=15,lwd=3,freq=FALSE,
         main=paste(obs$location,obs$obs.name),
         sub=paste(min(round(yy,2)),"--",max(round(yy,2)),
           ":",sub.tit,xlab=obs$unit))

    x.dist <- seq(min(histo$mids),max(histo$mids),length=101)
    y.dist <- dnorm(x.dist,
                       mean=mean(value,na.rm=TRUE),
                       sd=sd(value,na.rm=TRUE))
    
    lines(x.dist,y.dist,col="red")
    lines(x.dist,dgamma(x.dist-min(x.dist),
          shape=mean((value-min(x.dist))^2,na.rm=TRUE)/sd(value^2,na.rm=TRUE),
          scale=sd(value^2,na.rm=TRUE)/mean(value-min(x.dist),na.rm=TRUE)),
          col="blue",lty=3)
    grid()
    dev.off()
    file.copy(c(figname1,figname2),out.dir)
    file.remove(c(figname1,figname2))
  }

  plotStation <- list(yy=yy,mm=mm,value=value,loc=obs$location,
                      histo=histo,x.dist=x.dist,y.dist=y.dist)
  invisible(plotStation)
}
