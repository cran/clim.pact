# R.E. Benestad, met.no, Oslo, Norway 09.10.2002
# rasmus.benestad@met.no
#------------------------------------------------------------------------

newFig <- function() {
   dev <- paste(options()$device,"()",sep="")
   eval(parse(text=dev))
 }


plotDS <- function(ds.obj,leps=FALSE,plot.map=TRUE, plot.res=FALSE,
                   plot.rate=FALSE,add=FALSE,col="darkred",lwd=2,lty=1,
                   direc="output/") {

if (class(ds.obj)!="ds") stop("Need a 'ds' object!")
attach(ds.obj)

# Plotting: -----------------------------------------------
               
pred.descr <- paste("Empirical Downscaling (",id.1,"[")

lons <- lon.loc
lats <- lat.loc
for (i in 1:n.fld) {
  eval(parse(text=paste("x.srt<-order(lon.",i,")",sep="")))
  eval(parse(text=paste("y.srt<-order(lat.",i,")",sep="")))
  eval(parse(text=paste("ds.obj$lons<-ds.obj$lon.",i,"[x.srt]",sep="")))
  eval(parse(text=paste("ds.obj$lats<-ds.obj$lat.",i,"[y.srt]",sep="")))
  eval(parse(text=paste("ds.obj$X.",i,"<-ds.obj$X.",i,"[y.srt,x.srt]",sep="")))
  lons <- eval(parse(text=paste("c(lons,lon.",i,")",sep="")))
  lats <- eval(parse(text=paste("c(lats,lat.",i,")",sep="")))
}

subtitle <- paste("Calibration: ",month," ",v.name," at ",ds.obj$location,
                  " using ",id.1,": R2=",fit.r2,
                  "%, p-value=",fit.p,"%.",sep="")

print(paste("subtitle:",subtitle))


y.lim.tr <- range(c(y.o,pre.y,pre.gcm))
yymm.o<-yy.o + (mm.o-0.5)/12 + (dd.o-0.5)/365.25
yymm.gcm<-yy.gcm + (mm.gcm-0.5)/12 + (dd.gcm-0.5)/365.25

#if (!leps) par(ask=TRUE)
if ((!add) & (plot.map)) {
  if (leps) {
    figname<- paste("predictor_",v.name,"_",location,"_",region,"_",
                  month,ex.tag,".eps",sep="")
    postscript(file = figname,onefile=TRUE,horizontal=FALSE,paper="a4")
  } else eval(parse(text=paste(lower.case(options()$device),"()",sep="")))
  par(ps=16,cex.sub=0.7,cex.main=0.9)
  plot(c(floor(min(lons)),ceiling(max(lons))),
       c(floor(min(lats)),ceiling(max(lats))),type="n",
       main=paste(pred.descr,region,"] ->",v.name,")"),sub=subtitle,
       xlab="Time",ylab=paste(v.name,"(",unit,")"))

  col.tab=c("darkblue","darkred","darkgreen","brown")

  for (i in 1:n.fld) {
    eval(parse(text=paste("lines(c(min(lon.",i,"),max(lon.",i,")),",
                              "c(min(lat.",i,"),min(lat.",i,")),",
                              "col=col.tab[i],lty=2)",sep="")))
  eval(parse(text=paste("lines(c(min(lon.",i,"),max(lon.",i,")),",
                              "c(max(lat.",i,"),max(lat.",i,")),",
                              "col=col.tab[i],lty=2)",sep="")))
  eval(parse(text=paste("lines(c(min(lon.",i,"),min(lon.",i,")),",
                              "c(min(lat.",i,"),max(lat.",i,")),",
                              "col=col.tab[i],lty=2)",sep="")))
  eval(parse(text=paste("lines(c(max(lon.",i,"),max(lon.",i,")),",
                              "c(min(lat.",i,"),max(lat.",i,")),",
                              "col=col.tab[i],lty=2)",sep="")))
  eval(parse(text=paste("contour(lon.",i,",lat.",i,",t(ds.obj$X.",i,
               "),nlevels=7,add=TRUE,lwd=2,col=col.tab[i])",sep="")))
  }

  if (!is.null(lon.loc) & !is.null(lat.loc))
    points(lon.loc,lat.loc,pch=20,col="wheat",cex=1.5)
    points(lon.loc,lat.loc,pch=20,col="black",cex=0.9)
  addland()
  grid()

  if (n.fld > 1) {
    legend(min(c(lons,lon.loc)),
         max(c(lats,lat.loc)),
         c(pred.name[1:n.fld]),
         col=c(col.tab[1:n.fld]),
         lwd=2,lty=1,merge=TRUE,bg="grey95")
  }

  if (leps) { 
    dev.off()
    if (!file.exists(direc)){
      print(paste("The directory",direc,"does not exists.. Creates it.."))
      dir.create(direc)
    } 
    file.copy(figname,direc)
    file.remove(figname)
  } 
}

if (!add) {
  if (leps) {
    figname<- paste("scen_",v.name,"_",location,"_",region,"_",
                    month,ex.tag,".eps",sep="")
    postscript(file = figname,onefile=TRUE,horizontal=FALSE,paper="a4")
  } else newFig()
  par(ps=16,cex.sub=0.7,cex.main=0.9)

  plot(c(min(yymm.o[1],yymm.gcm[1]),yymm.gcm[length(yymm.gcm)]),
       y.lim.tr,type="n",
       main=paste(pred.descr,region,"] ->",v.name,")"),sub=subtitle,
       xlab="Time",ylab=paste(v.name,"(",unit,")"))
  grid()

}

lines(yymm.o,y.o,col="darkblue",lwd=3);
lines(yymm.o,pre.y,col="grey40",lty=2,lwd=2);
lines(yymm.gcm,pre.gcm,col=col,lwd=lwd,lty=lty);
lines(yymm.gcm,pre.fit, col = "red",lwd=1,lty=2) 
lines(yymm.gcm,pre.p.fit, col = "red",lwd=1,lty=2)
points(yymm.o,y.o,col="darkblue",pch=20);
points(yymm.o,pre.y,col="grey40",pch=21);
points(yymm.gcm,pre.gcm,col="darkred",pch=21);

if (!add) legend(quantile(c(yymm.o,yymm.gcm),0.01),
                 max(c(y.o,pre.y,pre.gcm)),
                 c("Obs.","Fit","GCM","Trends"),cex=0.75,
                 col=c("darkblue","grey40","darkred","red"),
                 lwd=c(3,2,2,1),lty=c(1,2,1,2),pch=c(20,21,21,26,26),
                 merge=TRUE,bg="grey95")

text(quantile(c(yymm.o,yymm.gcm),0.01),
     min(c(y.o,pre.y,pre.gcm)),pos=4,cex=0.6,
     paste(month,": Trend fit: P-value=",gcm.trnd.p,"%; ",
           "Projected trend= ",rate.ds,"+-",rate.err," ",
           unit,"/decade",sep=""))

if (leps) { 
  dev.off()
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
par(ps=16,cex.sub=0.7,cex.main=0.9)

plot(c(min(yymm.gcm),max(yymm.gcm)),y.lim.tr,type="n",
     main=paste(pred.descr,region,"] ->",
       v.name,")"),sub=subtitle,
     xlab="Time",ylab=paste("rate of change in",v.name,"(",unit,"/decade)"))
grid()
lines(yymm.gcm,tr.est.p.fit, col = "blue",lwd=3)
lines(c(min(yymm.gcm),max(yymm.gcm)),c(rate.ds,rate.ds),col = "red",lwd=2)

legend(min(yymm.gcm),-1.5,c("Polinomial fit","Linear fit"),
       lwd=c(3,2),col=c("blue","red"),bg="grey95")
if (leps) { 
  dev.off()
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
                region,"] ->",v.name,")"),sub=subtitle,
     xlab="Time",ylab=paste(v.name,"(",unit,")"))
lines(yymm.o,pre.y-mean(pre.y,na.rm=TRUE),col="grey",lty=3); 
grid()


if (leps) { 
  dev.off()
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
  dev.off()
  file.copy(figname,direc)
  file.remove(figname)
}
}
}
