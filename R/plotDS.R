# R.E. Benestad, met.no, Oslo, Norway 09.10.2002
# rasmus.benestad@met.no
#------------------------------------------------------------------------



plotDS <- function(ds.obj,leps=FALSE,plot.ts=TRUE,plot.map=TRUE, plot.res=FALSE,
                   plot.rate=FALSE,add=FALSE,col="darkred",lwd=2,lty=1,
                   direc="output/",main=NULL,sub=NULL,xlab=NULL,ylab=NULL) {

if (class(ds.obj)!="ds") stop("Need a 'ds' object!")
attach(ds.obj)

if (options()$device=="none") {
  plot.ts=FALSE; plot.map<- FALSE; plot.res<- FALSE; plot.rate <- FALSE
} 

# Plotting: -----------------------------------------------
  
pred.descr <- paste("Empirical Downscaling (",id.1,"[")

lons <- lon.loc
lats <- lat.loc
if (!is.finite(lons)) lons <- mean(ds.obj$lon,na.rm=TRUE)
if (!is.finite(lats)) lats <- mean(ds.obj$lat,na.rm=TRUE)

for (i in 1:n.fld) {
  eval(parse(text=paste("x.srt<-order(ds.obj$lon.",i,")",sep="")))
  eval(parse(text=paste("y.srt<-order(ds.obj$lat.",i,")",sep="")))
  eval(parse(text=paste("ds.obj$lons<-ds.obj$lon.",i,"[x.srt]",sep="")))
  eval(parse(text=paste("ds.obj$lats<-ds.obj$lat.",i,"[y.srt]",sep="")))
  ds.obj$X.1 <- ds.obj$X.1[y.srt,x.srt]
  eval(parse(text=paste("ds.obj$X.",i,"<-ds.obj$X.",i,"[y.srt,x.srt]",sep="")))
  lons <- eval(parse(text=paste("c(lons,lon.",i,")",sep="")))
  lats <- eval(parse(text=paste("c(lats,lat.",i,")",sep="")))
}

if (is.null(main)) main <-  paste(pred.descr,region,"] ->",v.name,")")       
if (is.null(sub)) sub <- paste("Calibration: ",month," ",v.name," at ",ds.obj$location,
                  " using ",id.1,": R2=",fit.r2,
                  "%, p-value=",fit.p,"%.",sep="")
if (is.null(xlab)) xlab <- "Time"
if (is.null(ylab)) ylab <- paste(v.name,"(",unit,")")

#print(paste("subtitle:",subtitle))

y.lim.tr <- range(c(ds.obj$y.o,ds.obj$pre.y,ds.obj$pre.gcm),na.rm=TRUE)
yymm.o<-ds.obj$yy.o + (ds.obj$mm.o-0.5)/12 + (ds.obj$dd.o-0.5)/365.25
yymm.gcm<-ds.obj$yy.gcm + (ds.obj$mm.gcm-0.5)/12 + (ds.obj$dd.gcm-0.5)/365.25

#if ((!leps)  & (plot.map)) par(ask=TRUE)
if (!add) {
  if ( (leps) & (plot.map) ) {
    figname<- paste("predictor_",v.name,"_",location,"_",region,"_",
                  month,ex.tag,".eps",sep="")
    postscript(file = figname,onefile=TRUE,horizontal=FALSE,paper="a4")
  } else newFig()
  if (plot.map) {
    par(ps=16,cex.sub=0.6,cex.main=0.7)
    plot(c(floor(min(lons,na.rm=TRUE)),ceiling(max(lons,na.rm=TRUE))),
         c(floor(min(lats,na.rm=TRUE)),ceiling(max(lats,na.rm=TRUE))),type="n",
         main=main,sub=sub,xlab=xlab,ylab=ylab)
  }

  col.tab=c("darkblue","darkred","darkgreen","brown")
  t.rng <- paste(range(ds.obj$yy.cal)[1],"-",range(ds.obj$yy.cal)[2])
  ds.map <- list(tim=NULL,date=NULL,n.maps=NULL)
  ds.map$tim<-month; ds.map$date<-t.rng; ds.map$n.maps=n.fld
  for (i in 1:n.fld) {
    if (plot.map) {
      eval(parse(text=paste("lines(c(min(ds.obj$lon.",i,"),max(ds.obj$lon.",i,")),",
                                "c(min(ds.obj$lat.",i,"),min(ds.obj$lat.",i,")),",
                                "col=col.tab[i],lty=2)",sep="")))
      eval(parse(text=paste("lines(c(min(ds.obj$lon.",i,"),max(ds.obj$lon.",i,")),",
                                "c(max(ds.obj$lat.",i,"),max(ds.obj$lat.",i,")),",
                                "col=col.tab[i],lty=2)",sep="")))
      eval(parse(text=paste("lines(c(min(ds.obj$lon.",i,"),min(ds.obj$lon.",i,")),",
                                "c(min(ds.obj$lat.",i,"),max(ds.obj$lat.",i,")),",
                                "col=col.tab[i],lty=2)",sep="")))
      eval(parse(text=paste("lines(c(max(ds.obj$lon.",i,"),max(ds.obj$lon.",i,")),",
                                "c(min(ds.obj$lat.",i,"),max(ds.obj$lat.",i,")),",
                                "col=col.tab[i],lty=2)",sep="")))
      eval(parse(text=paste("contour(ds.obj$lon.",i,",ds.obj$lat.",i,",t(ds.obj$X.",i,
                 "),nlevels=7,add=TRUE,lwd=2,col=col.tab[i])",sep="")))
    }
    eval(parse(text=paste("ds.map$lon.",i,"<-ds.obj$lon.",i,sep="")))
    eval(parse(text=paste("ds.map$lat.",i,"<-ds.obj$lat.",i,sep="")))
    eval(parse(text=paste("ds.map$map.",i,"<-t(ds.obj$X.",i,")",sep="")))
  }
  class(ds.map) <- "map"; attr(ds.map,"descr") <- "ds: large-scale pattern"

#  print("plotDS: HERE")

  if (plot.map) {
    if (!is.null(lon.loc) & !is.null(lat.loc))
      points(lon.loc,lat.loc,pch=20,col="wheat",cex=1.5)
      points(lon.loc,lat.loc,pch=20,col="black",cex=0.9)
    addland()
    grid()

    if (n.fld > 1) {
      legend(min(c(lons,lon.loc,na.rm=TRUE)),
           max(c(lats,lat.loc,na.rm=TRUE)),
           c(pred.name[1:n.fld]),
           col=c(col.tab[1:n.fld]),
           lwd=2,lty=1,merge=TRUE,bg="grey95")
    }

    if (leps) { 
        if (dev.cur()>1) dev.off()
      if (!file.exists(direc)){
        print(paste("The directory",direc,"does not exists.. Creates it.."))
        dir.create(direc)
      } 
      file.copy(figname,direc)
      file.remove(figname)
    } 
  }

}

if ((!add) & (plot.ts)) {
  if (leps) {
    figname<- paste("scen_",v.name,"_",location,"_",region,"_",
                    month,ex.tag,".eps",sep="")
    postscript(file = figname,onefile=TRUE,horizontal=FALSE,paper="a4")
  } else newFig()
  par(ps=16,cex.sub=0.6,cex.main=0.7,cex.main=0.6)

  plot(c(min(yymm.o[1],yymm.gcm[1]),yymm.gcm[length(yymm.gcm)]),
       y.lim.tr,type="n",
       main=main,sub=sub,xlab=xlab,ylab=ylab)
  grid()

}

if (plot.ts) {
  lines(yymm.o,ds.obj$y.o,col="darkblue",lwd=3);
  lines(yymm.o,ds.obj$pre.y,col="grey40",lty=2,lwd=2);
  lines(yymm.gcm,ds.obj$pre.gcm,col=col,lwd=lwd,lty=lty);
  lines(yymm.gcm,ds.obj$pre.fit, col = "red",lwd=1,lty=2) 
  lines(yymm.gcm,ds.obj$pre.p.fit, col = "red",lwd=1,lty=2)
  points(yymm.o,ds.obj$y.o,col="darkblue",pch=20);
  points(yymm.o,ds.obj$pre.y,col="grey40",pch=21);
  points(yymm.gcm,ds.obj$pre.gcm,col="darkred",pch=21);

  if (!add) legend(quantile(c(yymm.o,yymm.gcm),0.01),
                 max(c(ds.obj$y.o,ds.obj$pre.ds.obj$y,pre.gcm)),
                 c("Obs.","Fit","GCM","Trends"),cex=0.65,
                 col=c("darkblue","grey40","darkred","red"),
                 lwd=c(3,2,2,1),lty=c(1,2,1,2),pch=c(20,21,21,26,26),
                 merge=TRUE,bg="grey95")

  text(quantile(c(yymm.o,yymm.gcm),0.01),
     min(c(y.o,pre.y,pre.gcm)),pos=4,cex=0.6,
     paste(month,": Trend fit: P-value=",gcm.trnd.p,"%; ",
           "Projected trend= ",rate.ds,"+-",rate.err," ",
           unit,"/decade",sep=""))
}

if (leps) { 
  if (dev.cur()>1) dev.off()
  file.copy(figname,direc)
  file.remove(figname)
}

# Plot the rate of change:

if ((plot.rate) & !(add)) {
if (leps) { 
  figname<- paste("tendency_",v.name,"_",location,"_",region,"_",
                month,ex.tag,".eps",sep="")
  postscript(file = figname,onefile=TRUE,horizontal=FALSE,paper="a4")
} else newFig()
par(ps=16,cex.sub=0.6,cex.main=0.7)

plot(c(min(yymm.gcm),max(yymm.gcm)),y.lim.tr,type="n",
     main=main,sub=sub,xlab=xlab,ylab=paste("rate of change in",ylab))
grid()
lines(yymm.gcm,tr.est.p.fit, col = "blue",lwd=3)
lines(c(min(yymm.gcm),max(yymm.gcm)),c(rate.ds,rate.ds),col = "red",lwd=2)

legend(min(yymm.gcm),max(y.lim.tr),c("Polinomial fit","Linear fit"),
       lwd=c(3,2),col=c("blue","red"),cex=0.7,bg="grey95")
if (leps) { 
    if (dev.cur()>1) dev.off()
  file.copy(figname,direc)
  file.remove(figname)
}
}

# Plot the residuals:

if ((plot.res)  & !(add)) {
if (leps) { 
  figname<- paste("residual_",v.name,"_",location,"_",region,"_",
                month,ex.tag,".eps",sep="")
  postscript(file = figname,onefile=TRUE,horizontal=FALSE,paper="a4")
} else newFig()
par(ps=16,cex.sub=0.9,cex.main=0.7)
plot(yymm.o,step.wise$residual,type="l",lwd=3,
     main=paste("Residual",
                region,"] ->",v.name,")"),sub=sub,xlab=xlab,ylab=ylab)
lines(yymm.o,pre.y-mean(pre.y,na.rm=TRUE),col="grey",lty=3); 
grid()


if (leps) { 
    if (dev.cur()>1) dev.off()
  file.copy(figname,direc)
  file.remove(figname)
  figname<- paste("qq-residual_",v.name,"_",location,"_",region,"_",
                month,".eps",ex.tag,sep="")
  postscript(file = figname,onefile=TRUE,horizontal=FALSE,paper="a4")
} else newFig()
par(ps=16,cex.sub=0.8,cex.main=0.85)
qqnorm((step.wise$residual-mean(step.wise$residual,na.rm=TRUE))/
       sd(step.wise$residual,na.rm=TRUE))
lines(c(-5,5),c(-5,5),col="grey",lty=2)
grid()

if (leps) { 
    if (dev.cur()>1) dev.off()
  file.copy(figname,direc)
  file.remove(figname)
}
}

detach(ds.obj)
invisible(ds.map)
}
