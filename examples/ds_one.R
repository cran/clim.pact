rm(list=ls())
library(clim.pact)
source("KDVH.R")
source("clim.pact/R/plotDS.R"); source("clim.pact/R/plotStation.R")
source("clim.pact/R/plotEOF.R"); source("clim.pact/R/objDS.R")
source("clim.pact/R/corField.R"); source("clim.pact/R/eof.R")
source("clim.pact/R/ds.R"); source("clim.pact/R/avail.locs.R")
source("clim.pact/R/retrieve.nc.R"); source("clim.pact/R/cdfcont.R")
source("clim.pact/R/plotStation.R"); source("clim.pact/R/mergeStation.R")
source("clim.pact/R/anomaly.station.R"); source("clim.pact/R/getnarp.R")

ds.one <- function(ele=101,cmons=1:12,silent=TRUE,new.plot=TRUE,do.20c=TRUE,
                   do.a1b=TRUE,do.RCM=FALSE,qc=FALSE,scen="sresa1b",post=TRUE,
                   predictand = "narp",station="Tasiilaq",test=FALSE,off=FALSE) {

if (off) {
  ele <- 101; cmons <- 1:12; qc <- FALSE; silent <- TRUE; test <- TRUE
  new.plot <- TRUE; do.20c <- TRUE; do.a1b <- TRUE; do.rcm <- FALSE
  scen<- "sresa1b"; predictand <- "narp"; station <-"Tasiilaq"; post<-TRUE
}

v.nam <- switch(as.character(ele), '101'='tas', '601'='pr')
print("Get ERA40")
if (ele==101) {
  load("ERA40_t2m_mon.Rdata")
  era40.t2m <- catFields(t2m,lon=c(-90,50),lat=c(40,75)); rm(t2m)
} else if (ele==601) {
  load("ERA40_prec_mon.Rdata")
  era40.t2m <- catFields(prec,lon=c(-90,50),lat=c(40,75)); rm(prec)
  era40.t2m$dat <- era40.t2m$dat * 1.599876e-05/0.3675309
}

#print("Get NCEP")
#ncep.t2m <- retrieve.nc("/home/rasmusb/data/ncep/air.mon.mean.nc")

dir <- "/home/rasmusb/data/ipcc_FoAR/"
op.path <- "output/"
gcms <- list.files(path=dir,pattern=v.nam)
gcms <- gcms[grep(".nc",gcms)]
gcms <- gcms[grep(v.nam,gcms)]
gcms1 <- gcms[grep("20c3m",lower.case(gcms))]
gcms2 <- gcms[grep(scen,gcms)]

if (test) {gcms1 <- gcms1[1];gcms2 <- gcms2[1]}

print(gcms1)
print(gcms2)

#stations <- stations[9] # Tasilaq


if (lower.case(predictand)=="nordklim+metno") {
  obs1 <- getnordklim(station,ele=ele)
  param <- switch(as.character(ele),"101"="TAM","601"="RR")
  obs2 <- KDVH(obs1$station,param=param)
  obs <- mergeStation(obs1,obs2)
} else if (lower.case(predictand)=="nordklim") {
  obs <- getnordklim(station,ele=ele)
} else if (lower.case(predictand)=="nacd") {
  obs <- getnordklim(station,ele=ele) 
} else if (lower.case(predictand)=="narp") {
  obs <- getnarp(station,ele=ele)
  #print(summary(obs$val))
}

x.rng <- c(max(c(obs$lon-40,-180)),min(c(obs$lon+40,180)))
y.rng <- c(max(c(obs$lat-40,-90)), min(c(obs$lat+40,90)))
print(x.rng); print(y.rng)

plot(c(1890,2100),c(min(rowMeans(obs$val[,cmons]),na.rm=TRUE),
     2.5*max(rowMeans(obs$val[,cmons]),na.rm=TRUE)),
     type="n",main=obs$location,ylab=obs$obs.name,xlab="time")
grid()
obs.ts <- plotStation(obs,what="t",add=TRUE,col="grey20",type="p",pch=19,
                      l.anom=FALSE,mon=cmons,trend=TRUE,std.lev=FALSE)

i.gcm <- 0
for (gcm in gcms1) {
  i.gcm <- i.gcm + 1
  print(gcm)
  GCM <- retrieve.nc(paste(dir,gcm,sep=""),x.rng=x.rng,y.rng=y.rng,v.nam=v.nam,silent=TRUE)
  class(GCM) <- c("field","monthly.field.object")
  attr(GCM$tim,"unit") <- "month"; GCM$dd[] <- 15
  
  dot <- instring(".",gcm)
  gcm.nm <- substr(gcm,dot[2]+1,dot[2]+5)
  slsh <- instring("/",obs$location)
  if (slsh[1] > 0) {
    obs$location <- substr(obs$location,1,slsh[1]-1)
  }
  fname <- paste(op.path,"ds_one_4AR.",predictand,strip(obs$location),
                 obs$station,ele,".",gcm.nm,"c20",sep="")
  print(fname)

  ds.station <- objDS(era40.t2m,GCM,obs,plot=FALSE,qualitycontrol=qc,silent=silent)
  x <- ds.station$station
  plotStation(x,what="t",add=TRUE,col="grey40",type="l",lwd=2,
             lty=1,l.anom=FALSE,mon=cmons,trend=FALSE,std.lev=FALSE)
  text(x$yy[length(x$yy)],mean(x$val[length(x$yy),]),gcm.nm)
  plotStation(obs,what="t",add=TRUE,col="grey20",type="p",pch=19,
              l.anom=FALSE,mon=cmons,trend=TRUE,std.lev=FALSE)
  plotStation(obs,what="t",add=TRUE,col="grey20",type="l",lwd=1,lty=3,
              l.anom=FALSE,mon=cmons,trend=TRUE,std.lev=FALSE)
  print(paste("Saving in",fname))
  save(file=paste(fname,"Rdata",sep=""),x,ds.station)
  ds.scen <- data.frame(Year=x$yy,
                         Jan=round(x$val[,1],2),Feb=round(x$val[,2],2),
                         Mar=round(x$val[,3],2),Apr=round(x$val[,4],2),
                         May=round(x$val[,5],2),Jun=round(x$val[,6],2),
                         Jul=round(x$val[,7],2),Aug=round(x$val[,8],2),
                         Sep=round(x$val[,9],2),Oct=round(x$val[,10],2),
                         Nov=round(x$val[,11],2),Dec=round(x$val[,12],2))
  write.table(ds.scen,file=paste(fname,"txt",sep=""),row.names = FALSE,quote = FALSE, sep="\t ")
}


i.gcm <- 0
for (gcm in gcms2) {
  i.gcm <- i.gcm + 1
  print(gcm)
  GCM <- retrieve.nc(paste(dir,gcm,sep=""),x.rng=x.rng,y.rng=y.rng,v.nam=v.nam,silent=TRUE)
  class(GCM) <- c("field","monthly.field.object")
  attr(GCM$tim,"unit") <- "month"; GCM$dd[] <- 15
  
  dot <- instring(".",gcm)
  slsh <- instring("/",obs$location)
  if (slsh[1] > 0) {
    obs$location <- substr(obs$location,1,slsh[1]-1)
  }
  fname <- paste(op.path,"ds_4AR.",strip(obs$location),obs$station,
                 ".",substr(gcm,1,dot),"a1b",sep="")
  print(fname)
  ds.station <- objDS(era40.t2m,GCM,obs,plot=FALSE,qualitycontrol=qc,silent=silent)
print("HERE!")
  x <- ds.station$station
  x.ts <- plotStation(x,what="n",add=TRUE,col="steelblue",type="l",lwd=2,
                      lty=1,l.anom=FALSE,mon=cmons,trend=TRUE,std.lev=FALSE)
  x$val <- x$val - x.ts$trend[1] + obs.ts$trend[length(obs.ts$trend)] 
  plotStation(x,what="t",add=TRUE,col="steelblue",type="l",lwd=2,
              lty=1,l.anom=FALSE,mon=cmons,trend=TRUE,std.lev=FALSE)
  print(paste("Saving in",fname))
  save(file=paste(fname,"Rdata",sep=""),x,ds.station)
  ds.scen <- data.frame(Year=x$yy,
                         Jan=round(x$val[,1],2),Feb=round(x$val[,2],2),
                         Mar=round(x$val[,3],2),Apr=round(x$val[,4],2),
                         May=round(x$val[,5],2),Jun=round(x$val[,6],2),
                         Jul=round(x$val[,7],2),Aug=round(x$val[,8],2),
                         Sep=round(x$val[,9],2),Oct=round(x$val[,10],2),
                         Nov=round(x$val[,11],2),Dec=round(x$val[,12],2))
  write.table(ds.scen,file=paste(fname,"txt",sep=""),row.names = FALSE,quote = FALSE, sep="\t ")
}

plotStation(obs,what="t",add=TRUE,col="grey20",type="p",pch=19,
            l.anom=FALSE,mon=cmons,trend=TRUE,std.lev=FALSE)


if (do.rcm) {
 load("metno3_RCMs_oslo.Rdata")
 is <- 1
 yymm.RCM <- seq(2071+1/24,2100,length=360)
 ts.HCA2.am <- c(rep(NA,6),diff(cumsum(ts.HCA2[,is]),lag=12)/12,rep(NA,6))
 ts.HCB2.am <- c(rep(NA,6),diff(cumsum(ts.HCB2[,is]),lag=12)/12,rep(NA,6))
 ts.MPIA2.am <- c(rep(NA,6),diff(cumsum(ts.MPIA2[,is]),lag=12)/12,rep(NA,6))
 ts.MPIB2.am <- c(rep(NA,6),diff(cumsum(ts.MPIB2[,is]),lag=12)/12,rep(NA,6))
 points(yymm.RCM, ts.HCA2.am,pch=19,col="red",cex=0.8)
 points(yymm.RCM, ts.HCB2.am,pch=19,col="darkred",cex=0.8)
 points(yymm.RCM, ts.MPIA2.am,pch=19,col="blue",cex=0.8)
 points(yymm.RCM, ts.MPIB2.am,pch=19,col="darkblue",cex=0.8)
 lines(yymm.RCM, ts.HCA2.am,lty=3,col="red")
 lines(yymm.RCM, ts.HCB2.am,lty=3,col="darkred")
 lines(yymm.RCM, ts.MPIA2.am,lty=3,col="blue")
 lines(yymm.RCM, ts.MPIB2.am,lty=3,col="darkblue")

 legend(1895,2*max(rowMeans(obs$val[,cmons]),na.rm=TRUE),
        c("Obs.","E-DS 20C3M","E-DS A1b","RCM HCA2","RCM HCB2","RCM MPIA2","RCM MPIB2"),
        lty=c(0,1,1,0,0,0,0),lwd=c(0,1,1,1,1,1,1),pch=c(19,26,26,19,19,19,19),
        col=c("black","grey60","steelblue","red","darkred","blue","darkblue"),
        bg="grey97",cex=0.8)
}



  dev.copy2eps(file=paste("ds_one_",obs$location,obs$ele,".eps",sep=""))
  dev2bitmap(file=paste("ds_one_",obs$location,obs$ele,".png",sep=""))

 if (post) system(paste("/home/rasmusb/scripts/putftp ds_one_",
                  obs$location,obs$ele,".png",sep=""))
}

