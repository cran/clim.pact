plotDSobj <- function(result,outdir="output",figs=c(1,2,3,4),main="") {

if (class(result)!="objDS") {
  stop("The argument is not an 'objDS' object!")
}

months<-c("Jan","Feb","Mar","Apr","May","Jun",
          "Jul","Aug","Sep","Oct","Nov","Dec")

for (mon in months) {
   var.n <- paste("result$",mon,"$pre.gcm",sep="")
   eval(parse(text=paste(var.n,"<-",var.n," - mean(",var.n,",na.rm=TRUE)",sep="")))
   var.n <- paste("result$",mon,"$pre.y",sep="")
   eval(parse(text=paste(var.n,"<-",var.n," - mean(",var.n,",na.rm=TRUE)",sep="")))
   var.n <- paste("result$",mon,"$y.o",sep="")
   eval(parse(text=paste(var.n,"<-",var.n," - mean(",var.n,",na.rm=TRUE)",sep="")))
}

# Plotting and diagnostics:

# Construct a time series for the whole year:

  ds.all.gcm <-cbind(
           result$Jan$pre.gcm,result$Feb$pre.gcm,result$Mar$pre.gcm,
           result$Apr$pre.gcm,result$May$pre.gcm,result$Jun$pre.gcm,
           result$Jul$pre.gcm,result$Aug$pre.gcm,result$Sep$pre.gcm,
           result$Oct$pre.gcm,result$Nov$pre.gcm,result$Dec$pre.gcm)
  yymm.all.gcm <-cbind(
           result$Jan$yy.gcm + (result$Jan$mm.gcm - 0.5)/12,
           result$Feb$yy.gcm + (result$Feb$mm.gcm - 0.5)/12,
           result$Mar$yy.gcm + (result$Mar$mm.gcm - 0.5)/12,
           result$Apr$yy.gcm + (result$Apr$mm.gcm - 0.5)/12,
           result$May$yy.gcm + (result$May$mm.gcm - 0.5)/12,
           result$Jun$yy.gcm + (result$Jun$mm.gcm - 0.5)/12,
           result$Jul$yy.gcm + (result$Jul$mm.gcm - 0.5)/12,
           result$Aug$yy.gcm + (result$Aug$mm.gcm - 0.5)/12,
           result$Sep$yy.gcm + (result$Sep$mm.gcm - 0.5)/12,
           result$Oct$yy.gcm + (result$Oct$mm.gcm - 0.5)/12,
           result$Nov$yy.gcm + (result$Nov$mm.gcm - 0.5)/12,
           result$Dec$yy.gcm + (result$Dec$mm.gcm - 0.5)/12)
  y.gcm <- as.vector(t(ds.all.gcm))
  yymm.gcm <-  as.vector(t(yymm.all.gcm))
  ibad <- c(1,diff(yymm.gcm)) < 0
  yymm.gcm[ibad] <- NA

  ds.all.cal <-cbind(
           result$Jan$pre.y,result$Feb$pre.y,result$Mar$pre.y,
           result$Apr$pre.y,result$May$pre.y,result$Jun$pre.y,
           result$Jul$pre.y,result$Aug$pre.y,result$Sep$pre.y,
           result$Oct$pre.y,result$Nov$pre.y,result$Dec$pre.y)
  yymm.all.cal <-cbind(
           result$Jan$yy.cal + (result$Jan$mm.cal - 0.5)/12,
           result$Feb$yy.cal + (result$Feb$mm.cal - 0.5)/12,
           result$Mar$yy.cal + (result$Mar$mm.cal - 0.5)/12,
           result$Apr$yy.cal + (result$Apr$mm.cal - 0.5)/12,
           result$May$yy.cal + (result$May$mm.cal - 0.5)/12,
           result$Jun$yy.cal + (result$Jun$mm.cal - 0.5)/12,
           result$Jul$yy.cal + (result$Jul$mm.cal - 0.5)/12,
           result$Aug$yy.cal + (result$Aug$mm.cal - 0.5)/12,
           result$Sep$yy.cal + (result$Sep$mm.cal - 0.5)/12,
           result$Oct$yy.cal + (result$Oct$mm.cal - 0.5)/12,
           result$Nov$yy.cal + (result$Nov$mm.cal - 0.5)/12,
           result$Dec$yy.cal + (result$Dec$mm.cal - 0.5)/12)
  y.cal <- as.vector(t(ds.all.cal))
  yymm.cal <-  as.vector(t(yymm.all.cal))

  obs.all <-cbind(
           result$Jan$y.o,result$Feb$y.o,result$Mar$y.o,
           result$Apr$y.o,result$May$y.o,result$Jun$y.o,
           result$Jul$y.o,result$Aug$y.o,result$Sep$y.o,
           result$Oct$y.o,result$Nov$y.o,result$Dec$y.o)
  yymm.all.obs <-cbind(
           result$Jan$yy.o + (result$Jan$mm.o - 0.5)/12,
           result$Feb$yy.o + (result$Feb$mm.o - 0.5)/12,
           result$Mar$yy.o + (result$Mar$mm.o - 0.5)/12,
           result$Apr$yy.o + (result$Apr$mm.o - 0.5)/12,
           result$May$yy.o + (result$May$mm.o - 0.5)/12,
           result$Jun$yy.o + (result$Jun$mm.o - 0.5)/12,
           result$Jul$yy.o + (result$Jul$mm.o - 0.5)/12,
           result$Aug$yy.o + (result$Aug$mm.o - 0.5)/12,
           result$Sep$yy.o + (result$Sep$mm.o - 0.5)/12,
           result$Oct$yy.o + (result$Oct$mm.o - 0.5)/12,
           result$Nov$yy.o + (result$Nov$mm.o - 0.5)/12,
           result$Dec$yy.o + (result$Dec$mm.o - 0.5)/12)
  y.obs <- as.vector(t(obs.all))
  yymm.obs <-  as.vector(t(yymm.all.obs))

if (!is.null(result$Jan$f.name)) {
  slash<-instring("/",result$Jan$f.name)
  uscr<-instring("_",result$Jan$f.name)
  subtitle <- substr(result$Jan$f.name,slash+1,uscr[2]-1)
} else subtitle <- " "


if (sum(is.element(figs,1))>0) {newFig()
par(cex.main=0.8)

ele <- result$station$ele
if (ele==101) {
  val.rng <- c(-10,10)
  y.cal[y.cal > max(val.rng)] <- NA
  y.cal[y.cal < min(val.rng)] <- NA
  y.gcm[y.gcm > max(val.rng)] <- NA
  y.gcm[y.gcm < min(val.rng)] <- NA
} else if (ele==601) {
  val.rng <- c(0,300)
  y.cal[y.cal > max(val.rng)] <- NA
  y.cal[y.cal < min(val.rng)] <- NA
  y.gcm[y.gcm > max(val.rng)] <- NA
  y.gcm[y.gcm < min(val.rng)] <- NA
}                                
plot(range(yymm.obs,yymm.gcm,na.rm=TRUE),
     range(y.obs,y.gcm,y.cal,na.rm=TRUE),type="n",
     main=main,sub=subtitle,xlab="Time",ylab=result$Jan$unit)
#plot(range(yymm.obs,yymm.gcm,na.rm=TRUE),
#     range(y.obs,y.gcm,y.cal,na.rm=TRUE),type="n",
#     main=paste("Downscaled ",result$Jan$v.name," anomalies at ",result$Jan$location,
#                "     (",round(result$Jan$lat.loc,2),"N/",round(result$Jan$lon.loc,2),"E)",sep=""),
#     sub=subtitle,xlab="Time",ylab=result$Jan$unit)
grid()
points(yymm.obs+0.025,y.obs,pch=20,cex=1.2,col="grey60")
points(yymm.obs,y.obs,pch=20,cex=1.2,col="black")
lines(yymm.obs,y.obs,lty=3)
lines(yymm.cal+0.0025,y.cal,lty=2,lwd=2,col="grey60")
lines(yymm.cal,y.cal,lty=2,lwd=2,col="grey30")
lines(yymm.gcm+0.0025,y.gcm,lty=1,lwd=2,col="grey30")
lines(yymm.gcm,y.gcm,lty=1,lwd=1,col="grey50")
points(yymm.obs+0.025,y.obs,pch=20,cex=1.2,col="grey60")
points(yymm.obs,y.obs,pch=20,cex=1.2,col="black")
legend(min(yymm.obs,yymm.gcm,na.rm=TRUE),
       max(y.obs,y.gcm,y.cal,na.rm=TRUE),
       c("Observations                  ","Calibr.                ","Scenario                   "),
       col=c("black","grey30","grey50"),
       pch=c(20,26,26),lty=c(3,2,1),lwd=c(1,2,1),bg="grey95",cex=0.7)
if (lower.case(options()$device)=="x11") 
       dev.copy2eps(file=paste(outdir,"/plotDSobj_1.eps",sep=""))
}
# Residuals:

res.rng <- range(result$Jan$step.wise$residual,result$Feb$step.wise$residual,
                 result$Mar$step.wise$residual,result$Apr$step.wise$residual,
                 result$May$step.wise$residual,result$Jun$step.wise$residual,
                 result$Jul$step.wise$residual,result$Aug$step.wise$residual,
                 result$Sep$step.wise$residual,result$Oct$step.wise$residual,
                 result$Nov$step.wise$residual,result$Dec$step.wise$residual, na.rm=TRUE)
if (sum(is.element(figs,2))>0) {newFig()

par(cex.main=0.7)                                
plot(c(0,length(result$Jan$yy.o)),res.rng,     
     type="n", main=paste("Residuals ",result$Jan$v.name," anomalies at ",result$Jan$location,
                          "     (",round(result$Jan$lat.loc,2),"N/",round(result$Jan$lon.loc,2),"E)",sep=""),
     sub=subtitle,xlab="Time",ylab=result$Jan$unit)
grid()
points(result$Jan$step.wise$residual,col="black",pch=20,cex=0.5)
points(result$Feb$step.wise$residual,col="black",pch=20,cex=0.5)
points(result$Mar$step.wise$residual,col="black",pch=20,cex=0.5)
points(result$Apr$step.wise$residual,col="black",pch=20,cex=0.5)
points(result$May$step.wise$residual,col="black",pch=20,cex=0.5)
points(result$Jun$step.wise$residual,col="black",pch=20,cex=0.5)
points(result$Jul$step.wise$residual,col="black",pch=20,cex=0.5)
points(result$Aug$step.wise$residual,col="black",pch=20,cex=0.5)
points(result$Sep$step.wise$residual,col="black",pch=20,cex=0.5)
points(result$Oct$step.wise$residual,col="black",pch=20,cex=0.5)
points(result$Nov$step.wise$residual,col="black",pch=20,cex=0.5)
points(result$Dec$step.wise$residual,col="black",pch=20,cex=0.5)
lines(result$Jan$step.wise$residual,col="black")
lines(result$Feb$step.wise$residual,col="grey40")
lines(result$Mar$step.wise$residual,col="red")
lines(result$Apr$step.wise$residual,col="darkred")
lines(result$May$step.wise$residual,col="blue")
lines(result$Jun$step.wise$residual,col="darkblue")
lines(result$Jul$step.wise$residual,col="green")
lines(result$Aug$step.wise$residual,col="darkgreen")
lines(result$Sep$step.wise$residual,col="magenta")
lines(result$Oct$step.wise$residual,col="cyan")
lines(result$Nov$step.wise$residual,col="wheat")
lines(result$Dec$step.wise$residual,col="brown")
dev.copy2eps(file=paste(outdir,"/plotDSobj_2.eps",sep="")) 
}

if (sum(is.element(figs,3))>0) { newFig()
par(cex.main=0.8)
brks <- seq(res.rng[1]-1,res.rng[2]+1,length=25)
h.jan<-hist(result$Jan$step.wise$residual,breaks=brks)$density
h.feb<-hist(result$Feb$step.wise$residual,breaks=brks)$density
h.mar<-hist(result$Mar$step.wise$residual,breaks=brks)$density
h.apr<-hist(result$Apr$step.wise$residual,breaks=brks)$density
h.may<-hist(result$May$step.wise$residual,breaks=brks)$density
h.jun<-hist(result$Jun$step.wise$residual,breaks=brks)$density
h.jul<-hist(result$Jul$step.wise$residual,breaks=brks)$density
h.aug<-hist(result$Aug$step.wise$residual,breaks=brks)$density
h.sep<-hist(result$Sep$step.wise$residual,breaks=brks)$density
h.oct<-hist(result$Oct$step.wise$residual,breaks=brks)$density
h.nov<-hist(result$Nov$step.wise$residual,breaks=brks)$density
h.dec<-hist(result$Dec$step.wise$residual,breaks=brks)$density
brks <- hist(result$Dec$step.wise$residual,breaks=brks)$mids
plot(range(brks),range(c(h.jan,h.feb,h.mar,h.apr,h.may,h.jun,
     h.jul,h.aug,h.sep,h.oct,h.nov,h.dec)),type="n",
     main=paste("Residuals ",result$Jan$v.name," anomalies at ",result$Jan$location,
                          "     (",round(result$Jan$lat.loc,2),"N/",round(result$Jan$lon.loc,2),"E)",sep=""),
     sub=subtitle,ylab="Density",xlab=result$Jan$unit)
grid()
lines(brks,h.jan,col="black")
lines(brks,h.feb,col="grey40")
lines(brks,h.mar,col="red")
lines(brks,h.apr,col="darkred")
lines(brks,h.may,col="blue")
lines(brks,h.jun,col="darkblue")
lines(brks,h.jul,col="green")
lines(brks,h.aug,col="darkgreen")
lines(brks,h.sep,col="magenta")
lines(brks,h.oct,col="cyan")
lines(brks,h.nov,col="wheat")
lines(brks,h.dec,col="brown")
if (lower.case(options()$device)=="x11") 
     dev.copy2eps(file=paste(outdir,"/plotDSobj_3.eps",sep="")) 
}

rates <- c(result$Jan$rate.ds,result$Feb$rate.ds,result$Mar$rate.ds,
           result$Apr$rate.ds,result$May$rate.ds,result$Jun$rate.ds,
           result$Jul$rate.ds,result$Aug$rate.ds,result$Sep$rate.ds,
           result$Oct$rate.ds,result$Nov$rate.ds,result$Dec$rate.ds)
err <- c(result$Jan$rate.err,result$Feb$rate.err,result$Mar$rate.err,
           result$Apr$rate.err,result$May$rate.err,result$Jun$rate.err,
           result$Jul$rate.err,result$Aug$rate.err,result$Sep$rate.err,
           result$Oct$rate.err,result$Nov$rate.err,result$Dec$rate.err)
r2 <- c(result$Jan$fit.r2,result$Feb$fit.r2,result$Mar$fit.r2,
           result$Apr$fit.r2,result$May$fit.r2,result$Jun$fit.r2,
           result$Jul$fit.r2,result$Aug$fit.r2,result$Sep$fit.r2,
           result$Oct$fit.r2,result$Nov$fit.r2,result$Dec$fit.r2)
p.val <- as.numeric(c(result$Jan$gcm.trnd.p,result$Feb$gcm.trnd.p,result$Mar$gcm.trnd.p,
           result$Apr$gcm.trnd.p,result$May$gcm.trnd.p,result$Jun$gcm.trnd.p,
           result$Jul$gcm.trnd.p,result$Aug$gcm.trnd.p,result$Sep$gcm.trnd.p,
           result$Oct$gcm.trnd.p,result$Nov$gcm.trnd.p,result$Dec$gcm.trnd.p))

if (sum(is.element(figs,4))>0) {newFig()
par(col.axis="white",cex.main=0.8)
plot(c(0,25),range(c(rates+err,rates-err),na.rm=TRUE),type="n",
     main=paste("Linear trend rates ",result$Jan$v.name," derived ",result$Jan$location,
                          "     (",round(result$Jan$lat.loc,2),"N/",round(result$Jan$lon.loc,2),"E)",sep=""),
     sub=" ",ylab=paste(result$Jan$unit,"/ decade"),xlab="Month")
par(col.axis="black",ps=10,las=3)
axis(1, 1:24, rep(months,2))
axis(2)
grid()

scl <- diff(range(c(rates+err,rates-err),na.rm=TRUE))/10
polygon(c(1:24,reverse(1:24)),c(rep(rates+err,2),reverse(rep(rates-err,2))),
        col="wheat",border="grey",lwd=2)
lines(0:24+0.5,c(r2[1],rep(r2,2))/10*scl+min(rates-err,na.rm=TRUE),
      type="S",col="steelblue")
lines(rep(rates,2),lwd=2)
points((1:24)[rep(p.val < 5,2)],rep(rates[p.val < 5],2),pch=20,cex=1.5)
points((1:24)[rep(p.val >= 5,2)],rep(rates[p.val >= 5],2),pch=21,cex=1.5)
text((1:24)-0.33,rep(rates-0.01*diff(range(c(rates+err,rates-err),na.rm=TRUE)),2),rep(rates,2),
       cex=0.8,col="grey45")

for (i in 0:10) {
  lines(c(23.8,24),rep(i*scl + min(rates-err,na.rm=TRUE),2),col="steelblue")
  lines(c(0,24),rep(i*scl +    min(rates-err,na.rm=TRUE),2),lty=3,col="steelblue")
  text(23.5,i*scl +            min(rates-err,na.rm=TRUE),paste(i*10,'%',sep=""),
        cex=0.8,col="steelblue")
  }
mtext("R-squared (%) from calibration regression",side=4,col="steelblue",cex=0.80)
points(1,max(rates+err),pch=20); text(3,max(rates+err),"5% sign.level")
points(7,max(rates+err),pch=21); text(8,max(rates+err),"not sign.")

if (lower.case(options()$device)=="x11") 
      dev.copy2eps(file=paste(outdir,"/plotDSobj_4.eps",sep="")) 
}
}



objDS <- function(field.obs,field.gcm,station,plot=TRUE,positive=NULL,
                  mon=NULL,direc="output/",cal.id=NULL,
                  ldetrnd=TRUE,i.eofs=seq(1,8,by=1),ex.tag="",
                  method="lm",leps=FALSE,param="t2m",failure.action=NULL,
                  plot.res=FALSE,plot.rate=FALSE,xtr.args="",
                  swsm="step",predm="predict",lsave=FALSE,rmac=TRUE,
                  silent=FALSE,qualitycontrol=TRUE,LINPACK=TRUE,wOBS=0.25) {

  cmon<-c("Jan","Feb","Mar","Apr","May","Jun",
          "Jul","Aug","Sep","Oct","Nov","Dec")
  if (options()$device == "none") {plot <- FALSE; plot.res<- FALSE; plot.rate<- FALSE}

  dims <- dim(field.obs$dat); ny <- dims[2]; nx <- dims[3]
#  print(dims)
  wy <- 2*pi*seq(0,ny-1,by=1)/(ny-1)
  x.mod<-matrix(rep(NA,ny*nx),ny,nx)
  x.mod[,1]<-cos(wy);   x.mod[,2]<-sin(wy)
  x.mod[,3]<-cos(2*wy); x.mod[,4]<-sin(2*wy)
  x.mod[,5]<-cos(3*wy); x.mod[,6]<-sin(3*wy)
  x.mod[,7]<-cos(4*wy); x.mod[,8]<-sin(4*wy)
  wx <- 2*pi*seq(0,nx-1,by=1)/(nx-1)
  y.mod<-matrix(rep(NA,ny*nx),ny,nx)
  y.mod[1,]<-cos(wx);   y.mod[2,]<-sin(wx)
  y.mod[3,]<-cos(2*wx); y.mod[4,]<-sin(2*wx)
  y.mod[5,]<-cos(3*wx); y.mod[6,]<-sin(3*wx)
  y.mod[7,]<-cos(4*wx); y.mod[8,]<-sin(4*wx)
  if (is.null(mon)) mon  <-  1:12
  
  result <- list(station=station)
  print(paste("objDS: field.obs$v.name=",field.obs$v.name))
  if (is.null(positive) &
      sum(is.element(c("t2m","tem"),lower.case(substr(field.obs$v.name,1,3))))> 0) {
      positive <- TRUE
  }
  rates <- rep(NA,12)
  for (imon in mon) {
#    if (plot) {
#      newFig()
#      if (dev.cur() > 1) par(cex.sub=0.6,cex.axis=0.6,cex.lab=0.6,fin=c(2.37,2.37))
#    }
    #print(imon)
    cormap <- corField(field.obs,station,mon=imon,main="",plot=plot)
    if ((lower.case(options()$device)=="x11") & (plot)) {
       #print("HERE1"); print(dev.cur()); print(direc); print(options()$device)
       dev.copy2eps(file=paste(direc,"cormap_",cmon[imon],".eps",sep=""))
       #dev2bitmap(file=paste("cormap_",cmon[imon],".jpg",sep=""),type="jpeg",width=2.37,height=2.37,res=300)
       dev2bitmap(file=paste(direc,"cormap_",cmon[imon],".jpg",sep=""),res=300)
    }

    # Find optimal longitudes & latitudes:
    if (!silent) print(paste("No. points with correlation > 0.5=",sum(cormap$map>0.5,na.rm=TRUE)))
    if ( (sum(cormap$map>0.5,na.rm=TRUE)<1) & !is.null(failure.action) ) {
      if (!silent) print(">>> objDS: call screen.failure.action <<<")
      if (!silent) print(paste(":::  Maximum correlation=",max(cormap$map,na.rm=TRUE),
                                "    N. points r>0.5=",sum(cormap$map>0.5,na.rm=TRUE),
                                "    valid points=",sum(is.finite(cormap$map))))
      if (!silent) print(paste(failure.action,"(obs=",station,")",sep=""))
      ds <- eval(parse(text=paste(failure.action,"(obs=",station,", mon=",imon,")",sep="")))
    } else {
      
    latx <- 0.5*(field.obs$lat[2:ny]+field.obs$lat[1:(ny-1)])
    lonx <- 0.5*(field.obs$lon[2:nx]+field.obs$lon[1:(nx-1)])
    iy <- min( (1:ny)[station$lat <= field.obs$lat], na.rm=TRUE)
    ix <- min( (1:nx)[station$lon <= field.obs$lon], na.rm=TRUE)
    yprof <- as.vector(cormap$map[ix,]); yprof[is.na(yprof)] <- 0
    xprof <- as.vector(cormap$map[,iy]); xprof[is.na(xprof)] <- 0
    largescale <- data.frame(y=yprof, X=x.mod)      
    y.fit<-lm(y ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + X.7 + X.8,data=largescale)
    largescale <- data.frame(y=xprof, X=t(y.mod))
    lsX <- data.frame(X=x.mod);  lsY <- data.frame(X=t(y.mod));
    x.fit<-lm(y ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + X.7 + X.8,data=largescale)
    yhat <- predict(y.fit,newdata=lsX); xhat <- predict(x.fit,newdata=lsY);
    yzero <- yhat[2:ny]*yhat[1:(ny-1)]; xzero <- xhat[2:nx]*xhat[1:(nx-1)]
    lonx <- lonx[xzero < 0]; latx <- latx[yzero < 0]
    #print("TEST: new x.rng/y.rng:")
    x.rng <- c(max(c(min(field.obs$lon),max(lonx[lonx < station$lon])), na.rm=TRUE),
               min(c(max(field.obs$lon),min(lonx[lonx > station$lon])), na.rm=TRUE))
    y.rng <- c(max(c(min(field.obs$lat),max(latx[latx < station$lat])), na.rm=TRUE),
               min(c(max(field.obs$lat),min(latx[latx > station$lat])), na.rm=TRUE))
    #print(x.rng); print(y.rng)
    if (x.rng[1] > station$lon-10) x.rng[1] <- station$lon-10
    if (x.rng[2] < station$lon+10) x.rng[2] <- station$lon+10
    if (y.rng[1] > station$lat-10) y.rng[1] <- station$lat-10
    if (y.rng[2] < station$lat+10) y.rng[2] <- station$lat+10
    #print(x.rng); print(y.rng)
    if (plot) {
      #print("HERE1"); print(dev.cur()); print(direc); print(options()$device)
      par(cex.sub=0.6,cex.axis=0.6,cex.lab=0.6,fin=c(2.37,2.37),fig=c(0,1,0,1))
      plot(range(c(field.obs$lat,field.obs$lon)),range(cormap$map,na.rm=TRUE),type="n",
           main="",xlab="deg N & deg E",ylab="",
           sub=paste(round(field.obs$lon[ix],2),"E/",round(field.obs$lat[iy],2),"N",sep=""))
      grid()
      lines(range(c(field.obs$lat,field.obs$lon)),rep(0,2),lty=3)
      points(field.obs$lat,yprof)
      lines(field.obs$lat,as.numeric(yhat),lwd=2);
      lines(rep(y.rng[1],2),range(cormap$map,na.rm=TRUE),lty=2)
      lines(rep(y.rng[2],2),range(cormap$map,na.rm=TRUE),lty=2)

      points(field.obs$lon,xprof,col="red",pch=20)
      lines(field.obs$lon,as.numeric(xhat),col="red",lwd=2)
      lines(rep(x.rng[1],2),range(cormap$map,na.rm=TRUE),lty=2,col="red")
      lines(rep(x.rng[2],2),range(cormap$map,na.rm=TRUE),lty=2,col="red")
      if ((lower.case(options()$device)=="x11") & (plot)) {
        #print("HERE2"); print(dev.cur()); print(direc); print(options()$device)
        dev.copy2eps(file=paste(direc,"objDS_",cmon[imon],"_1.eps",sep=""))
        dev2bitmap(file=paste(direc,"objDS_",cmon[imon],"_1.jpg",sep=""),type="jpeg",width=2.37,height=2.37,res=300)
      }
    }
    print("catFields:")
    #print(">>> Check REB 11.02.2004!")
    #print(x.rng)
    #print(c(sum(!is.finite(field.obs$dat)),sum(!is.finite(field.gcm$dat))))
    #print(summary(field.obs$lon)); print(summary(field.obs$lat))
    #print(summary(field.gcm$lon)); print(summary(field.gcm$lat))
    #print(c(length(field.obs$yy),length(field.obs$mm),length(field.obs$id.t),NA,dim(field.obs$dat)))
    #print(c(length(field.gcm$yy),length(field.gcm$mm),length(field.gcm$id.t),NA,dim(field.gcm$dat)))

    field.2 <- catFields(field.obs,field.gcm,lon=x.rng,lat=y.rng,mon=imon)

    #print(c(length(field.2$yy),length(field.2$mm),length(field.2$id.t),NA,dim(field.2$dat)))
    if (!is.null(wOBS)) {                                    #REB 21.03.05
      if (!silent) print("Weight down GCM:")
      t.wgt <- wOBS*length(field.obs$id.t)/length(field.gcm$id.t)
      i.gcm <- is.element(field.2$id.t,field.gcm$id.t[1])
      field.2$dat[i.gcm,,] <- field.2$dat[i.gcm,,]*t.wgt
    }
    if (!silent) print("EOF:")
    eof <- EOF(field.2,silent=silent,plot=plot,lsave=FALSE,LINPACK=LINPACK)
    if (!is.null(wOBS)) {                                    #REB 21.03.05XS
      if (!silent) print("Weight up GCM:")
      i.gcm <- is.element(eof$id.t,field.gcm$id.t[1])
      eof$PC[i.gcm,] <- eof$PC[i.gcm,]/t.wgt
    }

#REB 09.03.05
#print(dim(eof$PC))
#print(length(eof$id.t))

    if (!silent) print("DS:")
    ds <- DS(preds=eof,dat=station,direc=direc,cal.id=cal.id,
                  ldetrnd=ldetrnd,i.eofs=i.eofs,ex.tag=ex.tag,
                  method=method,plot=FALSE,leps=leps,param=param,
                  plot.res=plot.res,plot.rate=plot.rate,xtr.args=xtr.args,
                  swsm=swsm,predm=predm,lsave=lsave,rmac=rmac,
                  silent=silent)
    if ( (ds$screening.failure) & !is.null(failure.action) ) {
      if (!silent) print(">>> objDS: call failure.action <<<")
      if (!silent) print(paste(failure.action,"(obs=station, mon=",imon,")",sep=""))
      ds <- eval(parse(text=paste(failure.action,"(obs=",station,", mon=",imon,")",sep="")))
    } else if (ds$screening.failure) print("screening failure - but continue as usual...")
  }
    ds$x.rng <- x.rng; ds$y.rng <- y.rng

    print("Grading for spatial pattern")
    field.x <- catFields(field.obs,lon=field.2$lon,lat=field.2$lat,mon=imon)
    dims <- dim(field.x$dat)
    if ( (length(dims)==3) & (dims[2]>1) & (dims[3]>1) ) {
      cormap <- corField(field.x,station,mon=imon,main="",plot=FALSE)
      patt.1 <- c(t(cormap$map))
      patt.2 <- c(ds$X.1) 
      valid <- is.finite(patt.1) & is.finite(patt.2)
    } else {
      print("objDS: Unexpected dimensions of field.x")
      print(dims)
      warning("objDS: Unexpected dimensions of field.x")
      patt.1 <- rep(0,length(ds$X.1))
      patt.2 <- c(ds$X.1) 
      valid <- is.finite(patt.2)      
    }
    if (sum(valid) > 30) grade.pattern <- round(10*cor(patt.1[valid],patt.2[valid]))/10 else
                         grade.pattern <- NA
    print(paste("Spatial correlation: ",sum(valid),"valid points. r=",grade.pattern))

    ds$grade.pattern <- grade.pattern

    #print(paste("result$",cmon[imon]," <- ds",sep=""))
    eval(parse(text=paste("result$",cmon[imon]," <- ds",sep="")))
    rates[imon] <- ds$rate.ds
  }
  
# Quality control!
# Check the rates with those of adjacent months: if very different, do the computation again
# but with a reduced domain size until a minimum size is reached (10x10 degrees).

#print(rates)
  drate <- -(diff(rep(rates,3))[11:22])*(diff(rep(rates,3))[12:23])
#print(drate)
  icheck <- (drate > 3*var(rates))
#print(icheck)  

  if (qualitycontrol) {
  print("==========================================================")
  print("=================== Quality control! =====================")
  print("==========================================================")
  while (sum(icheck)>0) {
    print(paste("Rates: (drate exceeds ",3*var(rates),")"))
    print(rbind(rates,drate,icheck))
    idoagain <- c((1:12)[icheck],(1:12)[icheck]+1,(1:12)[icheck]-1)
    idoagain[idoagain<1] <- 12; idoagain[idoagain > 12] <- 1
    for (ii in idoagain) {
       print(paste("Re-compute ",cmon[ii],": rate=",rates[ii]," drate=",drate[ii]))
       x.rng <- eval(parse(text=paste("result$",cmon[ii],"$x.rng",sep="")))
       y.rng <- eval(parse(text=paste("result$",cmon[ii],"$y.rng",sep="")))
       if (diff(x.rng) > 25) {
         if (station$lon > x.rng[1]+10) x.rng[1] <- x.rng[1]+(station$lon - x.rng[1] - 10)/3
         if (station$lon < x.rng[2]-10) x.rng[2] <- x.rng[2]-(x.rng[2] - station$lon - 10)/3
       }
       if (diff(y.rng) > 25) {
         if (station$lat > y.rng[1]+10) y.rng[1] <- y.rng[1]+(station$lat - y.rng[1] - 10)/3
         if (station$lat < y.rng[2]+10) y.rng[2] <- y.rng[2]-(y.rng[2] - station$lat - 10)/3
       }
       field.2 <- catFields(field.obs,field.gcm,lon=x.rng,lat=y.rng,mon=ii)
       if (!is.null(wOBS)) {                                    #REB 21.03.05
         if (!silent) print("Weight down GCM:")
         t.wgt <- wOBS*length(field.obs$id.t)/length(field.gcm$id.t)
         i.gcm <- is.element(field.2$id.t,field.gcm$id.t[1])
         field.2$dat[i.gcm,,] <- field.2$dat[i.gcm,,]*t.wgt
       }
       eof <- EOF(field.2,silent=TRUE,lsave=FALSE,LINPACK=LINPACK)
       if (!is.null(wOBS)) {                                    #REB 21.03.05XS
         if (!silent) print("Weight up GCM:")
         i.gcm <- is.element(eof$id.t,field.gcm$id.t[1])
         eof$PC[i.gcm,] <- eof$PC[i.gcm,]/t.wgt
       }
       ds <- DS(preds=eof,dat=station,direc=direc,cal.id=cal.id,
                ldetrnd=ldetrnd,i.eofs=i.eofs,ex.tag=ex.tag,
                method=method,plot=FALSE,leps=leps,param=param,
                plot.res=plot.res,plot.rate=plot.rate,xtr.args=xtr.args,
                swsm=swsm,predm=predm,lsave=lsave,rmac=rmac,
                silent=TRUE)
       ds$x.rng <- x.rng; ds$y.rng <- y.rng


       rates[ii] <- ds$rate.ds
       eval(parse(text=paste("result$",cmon[ii]," <- ds",sep="")))
       drate <- -(diff(rep(rates,3))[11:22])*(diff(rep(rates,3))[12:23])
       icheck[ii] <- (drate[ii] > 3*var(rates))
       if ((diff(x.rng) <= 15) & (diff(y.rng) <= 15)) icheck[ii] <- FALSE
       print(paste("New rate=",rates[ii]," drate=",drate[ii]," x.rng=",x.rng[1]," - ",x.rng[2],
                    " y.rng=",y.rng[1]," - ",y.rng[2],"icheck[ii]=",icheck[ii]))
    }
  }
  }

  grade.trend <- 1 - sum(icheck)/length(icheck)
  ds$grade.trend <- grade.trend

  if (!silent) print("The downscaling is complete for this station")
  class(result) <- "objDS"
 
  #print("ds.val:")
  ds.val <-cbind(
           result$Jan$pre.gcm,result$Feb$pre.gcm,result$Mar$pre.gcm,
           result$Apr$pre.gcm,result$May$pre.gcm,result$Jun$pre.gcm,
           result$Jul$pre.gcm,result$Aug$pre.gcm,result$Sep$pre.gcm,
           result$Oct$pre.gcm,result$Nov$pre.gcm,result$Dec$pre.gcm)
  #print("ds.yy:")
  
  ds.yy <- result$Jan$yy.gcm
  station.series <- station.obj(ds.val,ds.yy,obs.name=paste("downscaled",station$obs.name),
                                station$unit,ele=station$ele,mm=NULL,
                                station=station$station,lat=station$lat,
                                lon=station$lon,alt=station$alt,
                                location=station$location,wmo.no=station$wmo.no,
                                start=station$start,yy0=min(ds.yy),country=station$country,
   ref=paste("clim.pact >= v2.1-4 > objDS (Benestad, 2004, Eos, vol 85, #42, Oct 19, p.417):",
                                field.2$filename))
   result$station <- station.series

#print("objDS - HERE... plot?")
  if (plot) {
    plotDSobj(result,outdir=direc)
  }
#print("exit objDS")
  invisible(result)
}
