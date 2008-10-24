# R.E. Benestad, met.no, Oslo, Norway 22.05.2002
# rasmus.benestad@met.no
#-------------------------------------------------------------------
# PLot data from NORDKLIM stations.

plotStation <- function(obs,l.anom=TRUE,mon=NULL,
                        leps=FALSE,out.dir="output",what="b",trend=TRUE,std.lev=TRUE, 
                        type="b",pch=21,col="grey30",lwd=1,lty=1,add=FALSE,
                        main=NULL,sub=NULL,xlab=NULL,ylab=NULL,normal.period=NULL,method="rowMeans") {
#print(method)

if (sum(is.element(c("b","t","d","n"),what))==0) stop("Argumet 'what' must be 'b','t', 'n', or 'd'!")
if ( (class(obs)[2]!="monthly.station.record") &
     (class(obs)[2]!="daily.station.record") ){
  stop(paste("The predictand must be a 'monthly.station.record'",
             "object - Use  station.obj()"))
}

if (is.null(obs$ele)) {
  stop("The type of parameter cannot be determined: set 'ele' field (101 = temp., 601 = precip.)") 
}
if (is.null(main)) {
  if (class(obs)[2]=="monthly.station.record") main <- paste(obs$location,obs$obs.name) else 
                                               main <- obs$location
}
if (class(obs)[2]=="daily.station.record") {
  if (sum(is.finite(obs$t2m)) > 0) {
    if (!add) {
      newFig()
      plot(obs$yy + obs$mm/12 + obs$dd/365.25, obs$t2m,pch=20,cex=0.5,
           main=main,sub="met.no Klima DataVareHus",
           xlab="Time",ylab="Temperature (deg C)")
      grid()
    }
    lines(obs$yy + obs$mm/12 + obs$dd/365.25,obs$t2m,lty=3,col=col,type="b",pch=pch,cex=0.5)
  }

  if (sum(is.finite(obs$precip)) > 0) {
    newFig()
    plot(obs$yy + obs$mm/12 + obs$dd/365.25, obs$precip,pch=20,cex=0.5,
         main=main,sub="met.no Klima DataVareHus",
         xlab="Time",ylab="Precipitation (mm)")
    grid()
    lines(obs$yy + obs$mm/12 + obs$dd/365.25, obs$precip,lty=3,col="grey")
  }
  plotStation <- obs


} else if (class(obs)[2]=="monthly.station.record") {
  
if ((!obs$found) | (sum(is.finite(obs$val))==0)) stop("No valid data!")

cmon<-c("Jan","Feb","Mar","Apr","May","Jun",
        "Jul","Aug","Sep","Oct","Nov","Dec")

if (!is.null(mon)) {
  if (((length(mon)== 1)) & (mon[1] > 0)) season <- cmon[mon]
  if (((length(mon)== 1)) & (mon[1] ==0)) season <- ""
  if ((length(mon)> 1)) season <- paste(cmon[mon[1]],'-',
                                      cmon[mon[length(mon)]],sep="")
} else {
  season <-"Dec-Jan"
}

  loc <- obs$location
  ny <- length(obs$yy)
  mm <- rep(1:12,ny)

  obsa <- anomaly.station(obs,period=normal.period)

  if (is.null(mon)) {

    if (l.anom) {
      value <- t(obsa$val)
        
      if ((!is.null(obs$alt)) & (!is.null(obs$lon)) & (!is.null(obs$lat))) {
          sub.tit <- paste("Anomaly",round(obs$alt,2),"m a.sl.",
                         round(obs$lon,2),"degE",
                         round(obs$lat,2),"degN")
      } else sub.tit <- paste("Anomaly:",loc)
          
    } else {
     value <- t(obs$val)
     if ((!is.null(obs$alt)) & (!is.null(obs$lon)) & (!is.null(obs$lat))) {
         sub.tit <- paste("Absolute",round(obs$alt,2),"m a.sl.",
                          round(obs$lon,2),"degE",
                          round(obs$lat,2),"degN")
       } else sub.tit <- paste("Absolute:",loc)
    }
    dims <- dim(value)
    dim(value) <- c(dims[1]*dims[2],1)
    yymm <- sort(rep(obs$yy,12)) + (rep(seq(1,12,by=1),ny)-0.5)/12
    yy <- sort(rep(obs$yy,12)); mm <- rep(seq(1,12,by=1),ny)
    clim <- mean(obs$val[1,] - obsa$val[1,],na.rm=TRUE)

  } else {
    if (l.anom) obs <- obsa
    yy <- obs$yy
    ny <- length(obs$yy)
    mm <- rep(as.vector(mon[1]),ny)
    yymm <- yy + (mm -0.5)/12
    #print(dim(obs$val))
    #print(length(rowMeans(obs$val[,mon])))
    #print(length(colMeans(obs$val[,mon])))
    if ( (length(mon)==3) & (mon[1]==12) & (mon[2]==1) & (mon[3]==2)) {
      obs$val[2:ny,12] <- obs$val[1:(ny-1),12]
      obs$val[1,12] <- NA
    } else if ( (length(mon)==3) & (mon[1]==11) & (mon[2]==12) & (mon[3]==1)) {  # REB Fix 21.10.2005
      obs$val[2:ny,11:12] <- obs$val[1:(ny-1),11:12]
      obs$val[1,11:12] <- NA
    } 

    # REB 10.04.2008
    if (length(mon)>1) value <- eval(parse(text=paste(method,"(obs$val[,mon])",sep=""))) else
                       value <- obs$val[,mon]
#    if (is.element(obs$ele,c(101,111,121,401,601,701,801,911)))
#          for (i in 1:ny) value[i] <- mean(obs$val[i,mon],na.rm=TRUE)
    if (is.element(obs$ele,c(112,602)))
          for (i in 1:ny) value[i] <- max(obs$val[i,mon],na.rm=TRUE)
    if (is.element(obs$ele,c(122)))
          for (i in 1:ny) value[i] <- max(obs$val[i,mon],na.rm=TRUE)
    
    if ((!is.null(obs$alt)) & (!is.null(obs$lon)) & (!is.null(obs$lat))) {
       sub.tit <- paste(season," - ",round(obs$alt,2),"m a.sl.",
                     round(obs$lon,2),"degE",
                     round(obs$lat,2),"degN")
    } else sub.tit <- paste(season,": ",loc,sep="")
    clim <- mean(obs$val[1,mon] - obsa$val[1,mon],na.rm=TRUE)

  }
  stdv <- sd(value,na.rm=TRUE)

  # Polinomial trend

  if (trend) {
    y <- value
    x <- yymm
    X <- data.frame(y=value,x = yymm)
    lm.tr.p<-lm(y ~ x + I(x^2) +I(x^3) + I(x^4) + I(x^5))
    pre.p.fit<-predict(lm.tr.p,newdata=X)
    coef.p.fit<-lm.tr.p$coefficients
    coef.p.fit[is.na(coef.p.fit)] <- 0
    der.p.fit<-c(coef.p.fit[2],2*coef.p.fit[3],3*coef.p.fit[4],
               4*coef.p.fit[5],5*coef.p.fit[6])
    tr.est.p.fit<-(der.p.fit[1] + der.p.fit[2]*yymm + der.p.fit[3]*yymm^2 +
                   der.p.fit[4]*yymm^3 + der.p.fit[5]*yymm^4)*10
  } else pre.p.fit <- NULL

  good <- !is.na(value)
  YYMM <- yymm; YY <- yy; MM <- mm; VALUE <- value; PRE.p.fit <- pre.p.fit
  yymm <- yymm[good]; yy <- yy[good]; mm <- mm[good]
  pre.p.fit <- pre.p.fit[good]
  value <- value[good]

  if (!leps) {
    
#  par(ask=TRUE)
    if ((what=="t") | (what=="b")) {
      if (!add) {
        newFig()
        par(cex.sub=0.8)
        plot(yymm,value,type=type,lwd=lwd,col=col,pch=pch,lty=lty,
                     main=main,sub=sub.tit,xlab="Time",ylab=obs$unit,cex=0.6)
      } else {
        if (type!="p") lines(yymm,value,lwd=lwd,col=col,lty=lty) 
        if ((type=="p") | (type=="b")) points(yymm,value,col=col,pch=pch,cex=0.6) 
      }
      if (trend) lines(yymm,pre.p.fit,col="red") 
      if (std.lev) lines(c(min(yymm),max(yymm)),rep(mean(value,na.rm=TRUE)+
                                   1.96*sd(value,na.rm=TRUE),2),
                                   lty=2,col="grey")
      if (std.lev) lines(c(min(yymm),max(yymm)),rep(mean(value,na.rm=TRUE)-
                                   1.96*sd(value,na.rm=TRUE),2),
                                   lty=2,col="grey")
      grid()
    }

    if ((what=="d") | (what=="b")) {
      newFig()
      par(cex.sub=0.8)
      if (!is.finite(clim)) clim <- 0
      histo <- hist(value[!is.na(value)],breaks=15,lwd=3,
         main=main,
         sub=paste(min(round(yy,2)),"--",max(round(yy,2)),
           ":",sub.tit,xlab=obs$unit),xlab=paste(obs$obs.name,obs$unit))

      x.dist <- seq(min(histo$mids),max(histo$mids),length=101)
      y.dist <- dnorm(x.dist,
                      mean=mean(value,na.rm=TRUE),
                      sd=sd(value,na.rm=TRUE))
      if (trend) lines(x.dist,y.dist,col="red")
#      lines(x.dist,dgamma(x.dist-min(x.dist),
#            shape=mean((value-min(x.dist))^2,na.rm=TRUE)/sd(value^2,na.rm=TRUE),
#            scale=sd(value^2,na.rm=TRUE)/mean(value-min(x.dist),na.rm=TRUE)),
#            col="blue",lty=3)
      if (l.anom) {
        for (i in -5:5) {
          lines(rep(0+i*stdv,2),c(0,0.005),lwd=3,col="grey40")
          text(i*stdv,0.015,round(clim+i*stdv,1),font=2,col="grey40") 
        }
         lines(rep(0,2),c(0.025,0.05),lwd=3,col="grey40")
      } else {
        for (i in -5:5) {
          lines(rep(clim+i*stdv,2),c(0,0.005),lwd=3,col="grey40")
          text(clim+i*stdv,0.015,round(i*stdv,1),font=2,col="grey40")
      }
        lines(rep(clim,2),c(0.025,0.05),lwd=3,col="grey40")
      }
      grid()
    } else  {
      histo <- hist(value[!is.na(value)],breaks=15,lwd=3,plot=FALSE)
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
    plot(yymm,value,type=type,lwd=3,
         main=main,
         sub=sub.tit,xlab="Time",ylab=obs$unit)
    if (trend) lines(yymm[!is.na(y)],pre.p.fit,col="red")
    lines(c(min(yymm),max(yymm)),rep(mean(value,na.rm=TRUE)+
                                 1.96*sd(value,na.rm=TRUE),2),
          lty=2,col="grey")
    lines(c(min(yymm),max(yymm)),rep(mean(value,na.rm=TRUE)-
                                 1.96*sd(value,na.rm=TRUE),2),
          lty=2,col="grey")
    grid()
    dev.off()

    postscript(file = figname2,onefile=TRUE,horizontal=FALSE,paper="a4")
    par(ps=14,cex.sub=0.8)
    histo <- hist(value,breaks=15,lwd=3,
         main=main,
         sub=paste(min(round(yy,2)),"--",max(round(yy,2)),
           ":",sub.tit,xlab=obs$unit))

    x.dist <- seq(min(histo$mids),max(histo$mids),length=101)
    y.dist <- dnorm(x.dist,
                       mean=mean(value,na.rm=TRUE),
                       sd=sd(value,na.rm=TRUE))
    
    if (trend) lines(x.dist,y.dist,col="red")
    lines(x.dist,dgamma(x.dist-min(x.dist),
          shape=mean((value-min(x.dist))^2,na.rm=TRUE)/sd(value^2,na.rm=TRUE),
          scale=sd(value^2,na.rm=TRUE)/mean(value-min(x.dist),na.rm=TRUE)),
          col="blue",lty=3)
    grid()
    dev.off()
    file.copy(c(figname1,figname2),out.dir)
    file.remove(c(figname1,figname2))
  }

  plotStation <- list(yy=yy,mm=mm,yymm=yymm,value=value,loc=obs$location,
                      histo=histo,x.dist=x.dist,y.dist=y.dist,trend=PRE.p.fit)

}
  invisible(plotStation)

}
