# R.E. Benestad, met.no, Oslo, Norway 09.10.2002
# rasmus.benestad@met.no
#------------------------------------------------------------------------

plumePlot <- function(ds.name.list=NULL,location,mon,direc="output",
                         t.rng=c(1850,2074),r2.th=50,p.th=0.05,
                         col="darkred",lwd=2,lty=1) {

  cmon<-c('Jan','Feb','Mar','Apr','May','Jun',
        'Jul','Aug','Sep','Oct','Nov','Dec')
  
  if (is.null(ds.name.list)) ds.name.list <- avail.ds(direc=direc)

  yy <- seq(t.rng[1],t.rng[2],by=1)
  nt <- length(yy)
  nds <- length(ds.name.list)
  sce <- matrix(rep(NA,nt*nds),nt,nds)
  i.sce <- rep(FALSE,nds)
  yy.min <- NA
  yy.max <- NA
  
#  x11()
#  load(paste(direc,"/",ds.name.list[1],sep=""))
#  plot(ds$yy.gcm,ds$pre.gcm,type="n")
  
  for (i.ds in 1:nds) {
    load(paste(direc,"/",ds.name.list[i.ds],sep=""))
#    print(summary(ds))
#    print(c(strip(lower.case(ds$location)),lower.case(location)))
#    print(ds$mm.gcm[1])
#    print(c(ds$fit.r2,ds$fit.p))
    if ( (strip(lower.case(ds$location))==lower.case(location)) &
         (mon==ds$mm.gcm[1]) & (ds$fit.r2 >= r2.th) &
         (ds$fit.p <= p.th) ) {
      i1 <- is.element(yy,ds$yy.gcm)
      i2 <- is.element(ds$yy.gcm,yy)
      sce[i1,i.ds] <- ds$pre.gcm[i2]
      i.sce[i.ds] <- TRUE
#      print(range(ds$yy.gcm[i2]))
#      lines(ds$yy.gcm[ii],ds$pre.gcm[ii])
      if ((is.na(yy.min)) | (min(ds$yy.gcm)<yy.min)) yy.min<-min(ds$yy.gcm) 
      if ((is.na(yy.max)) | (max(ds$yy.gcm)>yy.max)) yy.max<-max(ds$yy.gcm) 
    }
  }
  sce <- sce[,i.sce]
  q975 <- rep(NA,nt)
  q025 <- rep(NA,nt)
  q500 <- rep(NA,nt)
  q250 <- rep(NA,nt)
  q750 <- rep(NA,nt)
  for (it in 1:nt) {
    if (sum(is.finite(sce[it,]))>0) {
      q975[it] <- quantile(sce[it,is.finite(sce[it,])],0.975)
      q025[it] <- quantile(sce[it,is.finite(sce[it,])],0.025)
      q500[it] <- quantile(sce[it,is.finite(sce[it,])],0.500)
      q250[it] <- quantile(sce[it,is.finite(sce[it,])],0.250)
      q750[it] <- quantile(sce[it,is.finite(sce[it,])],0.750)
    }
  }
#  print(summary(sce))
#  print(i.sce)
  x.ind <- seq(0,1,length=nt)
  tr.dat<-data.frame(y=q975, x=x.ind)
  lm.tr.p<-lm(y ~ x + I(x^2) +I(x^3), data=tr.dat)
  p975<-predict(lm.tr.p,newdata=tr.dat)
  tr.dat<-data.frame(y=q025, x=x.ind)
  lm.tr.p<-lm(y ~ x + I(x^2) +I(x^3), data=tr.dat)
  p025<-predict(lm.tr.p,newdata=tr.dat)
  tr.dat<-data.frame(y=q500, x=x.ind)
  lm.tr.p<-lm(y ~ x + I(x^2) +I(x^3), data=tr.dat)
  p500<-predict(lm.tr.p,newdata=tr.dat)
  tr.dat<-data.frame(y=q250, x=x.ind)
  lm.tr.p<-lm(y ~ x + I(x^2) +I(x^3), data=tr.dat)
  p250<-predict(lm.tr.p,newdata=tr.dat)
  tr.dat<-data.frame(y=q750, x=x.ind)
  lm.tr.p<-lm(y ~ x + I(x^2) +I(x^3), data=tr.dat)
  p750<-predict(lm.tr.p,newdata=tr.dat)       

  subtitle <- paste("Based on",sum(i.sce),"scenarios.",
                    "Using a cubic-fit to trend,",
                    "R^2 fit >",r2.th,"% & p-val fit <",p.th)
  yy.o <- ds$yy.o
  obs <- ds$y.o
  ii <- is.element(yy,seq(yy.min,yy.max,by=1))
  y.lim.tr <- range(c(p975[ii],p025[ii]))

  print(summary(p025))
  plot(c(min(c(yy.o,yy.min)),yy.max),
       y.lim.tr,type="n",
       main=location,sub=subtitle,
       xlab="Time",
       ylab=paste(cmon[mon],ds$v.name,"(",ds$unit,")"))
  grid()

  polygon(c(yy,reverse(yy)),
          c(p975,reverse(p025)),col="grey75")
  polygon(c(yy,reverse(yy)),
          c(p250,reverse(p750)),col="wheat",density=10,lwd=5)
  lines(yy,p500,col="black",lwd=3)
  lines(yy,q975,col="grey20",type="s",lty=3)
  lines(yy,q025,col="grey20",type="s",lty=3)
  lines(yy,q250,col="grey20",type="s",lty=3)
  lines(yy,q750,col="grey20",type="s",lty=3)
  lines(yy,q500,col="black",type="s",lty=3)
  lines(range(yy),rep(min(p500[ii]),2),lty=2)
  lines(range(yy),rep(max(p500[ii]),2),lty=2)
  points(yy.o,obs,pch=20)
  points(yy.o,obs,pch=20,cex=0.6,col="red")

  legend(quantile(yy,0.95),
         quantile(c(p025[ii]),0.75),
         c("2.5%--97.5%","25%--75%","Median","Obs"),
         col=c("grey75","wheat","black","red"),
         lwd=c(5,5,3,0),lty=c(rep(1,3),0),
         pch=c(rep(26,3),20),bg="grey95",cex=0.6)
}
