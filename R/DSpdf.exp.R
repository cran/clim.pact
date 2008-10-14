
ExtEOF <- function(fields,lag=1,mon=NULL,lon=NULL,lat=NULL) {
  X.e <- fields
  id.name <- paste(X.e$id.t[1],0,sep="")
  X.e$id.x[,]<-id.name; X.e$id.x[]<-id.name; X.e$id.lon[]<-id.name; X.e$id.lat[]<-id.name
  if (lag==0) eeof <- EOF(fields,mon=mon,lon=lon,lat=lat) else {
    for (i in 1:length(lag)) {
      X.l <- lagField(fields,lag[i])
      id.name <- paste(X.e$id.t[1],i,sep="")
      X.l$id.x[,]<-id.name; X.l$id.x[]<-id.name; X.l$id.lon[]<-id.name; X.l$id.lat[]<-id.name
      X.e <-  mixFields(X.e,X.l)
    }
    eeof <- EOF(X.e,mon=mon,lon=lon,lat=lat)
  }
  invisible(eeof)
}


lagField <- function(fields,lag=1) {
  fields$mm <- fields$mm - lag
  zeros <- fields$mm==0
  if (sum(zeros)>0) {
    fields$mm[zeros] <- 12
    fields$yy[zeros] <- fields$yy[zeros] - 1
  }
  thirteens <- fields$mm==13
  if (sum(thirteens)>0) {
    fields$mm[thirteens] <- 1
    fields$yy[thirteens] <- fields$yy[zeros] + 1
  }
  invisible(fields)
}

DSpdf.exp <- function(obs=NULL,dT=0,dP=0,plot=TRUE,year=NULL,month=NULL) {
  data(exp.law1)
  data(addland)
  dist <- min(distAB(obs$lon,obs$lat,lon.cont,lat.cont),na.rm=TRUE)/1000
  
  slope <- data.frame(slope=exp.par$slope,temp=exp.par$mt2m,lon=exp.par$lons,
                      lat=exp.par$lats,alt=exp.par$alt,dist=exp.par$dist,
                      precip=exp.par$mprecip)
  const <- data.frame(const=exp.par$const,temp=exp.par$mt2m,lon=exp.par$lons,
                      lat=exp.par$lats,alt=exp.par$alt,dist=exp.par$dist,
                      precip=exp.par$mprecip)
  ii <- is.finite(obs$t2m) & is.finite(obs$precip)
  if (!is.null(year))  ii <- ii & is.element(obs$yy,year)
  if (!is.null(month)) ii <- ii & is.element(obs$mm,month)

  extrap.dep <- data.frame(temp=mean(obs$t2m[ii]),lon=obs$lon,lat=obs$lat,
                       alt=obs$alt,dist=dist,precip=mean(obs$precip[ii]))
  extrap.chg <- data.frame(temp=mean(obs$t2m[ii])+dT,lon=obs$lon,lat=obs$lat,
                         alt=obs$alt,dist=dist,precip=mean(obs$precip[ii])+dP)
  slope.model <- lm(slope ~ temp + precip + lon + lat + alt + dist,data=slope)
  const.model <- lm(const ~ temp + precip + lon + lat + alt + dist,data=const)
  smod <- step(slope.model,trace=0)
  cmod <- step(const.model,trace=0)
  slope.dep <- predict(smod,newdata=extrap.dep)
  const.dep <- predict(cmod,newdata=extrap.dep)
  slope.chg <- predict(smod,newdata=extrap.chg)
  const.chg <- predict(cmod,newdata=extrap.chg)
  
  x <- obs$precip[ii]; x <- x[x > min(exp.par$minAmountPrecip)]
  exp.y <- as.numeric(table(round(x)))
  exp.x <- as.numeric(rownames(table(round(x))))
  h <- exp.y /(sum(exp.y)*min(diff(exp.x))) 
  log.data <- data.frame(y=log(exp.y),x=exp.x)

  pdf <- exp(slope.dep*exp.x)/(sum(exp(slope.dep*exp.x))*min(diff(exp.x)))
  pdf.chg <- exp(slope.chg*exp.x)/(sum(exp(slope.chg*exp.x))*min(diff(exp.x)))
  model <- paste("exp[ ",round(slope.dep,4),"x ]")
  model.chg <- paste("exp[ ",round(slope.chg,4),"x ]")
  log.mod <- lm(y ~ x,data=log.data)
  log.fit <- const.dep + slope.dep*exp.x
  log.fit.chg <- const.chg + slope.chg*exp.x
  log.fit[log.fit < 0] <- NA; log.fit.chg[log.fit.chg < 0] <- NA

  if (plot) {
    plot(c(0,100),c(0,0.2),type="n",main=obs$location,sub="",
         xlab="Precipitation (mm/day)",ylab="density")
    grid()
    points(exp.x,h,pch=20,col="grey70",cex=1.5)
    lines(exp.x,pdf,lwd=2)

    polygon(c(50,100,100,50,50)+1,c(0.1,0.1,0.2,0.2,0.1)+0.002,col="grey70",
            border="grey90",lwd=2)
    polygon(c(50,100,100,50,50),c(0.1,0.1,0.2,0.2,0.1),col="grey97")
    points(exp.x/2+50,log(exp.y)/100+0.1,col="darkred",cex=0.9,pch=20)
    points(exp.x/2+50,log(exp.y)/100+0.1,col="red",cex=0.8,pch=20)
    lines(exp.x/2+50,log.fit/100+0.1,lty=2,lwd=2)

    text(48,0.15,"ln(density)",srt=90,cex=0.8)    
    for (i in seq(1,100,by=10)) lines(rep(i,2)/2+50,c(0.10,0.101))
    text(70,0.09,model,cex=0.8)
    text(70,0.07,paste("Low precip cut-off:"=exp.par$minAmountPrecip),cex=0.8)
    if (dT != 0) {
      lines(exp.x,pdf.chg,lty=2,lwd=1,col="steelblue")
      lines(exp.x/2+50,log.fit.chg/100+0.1,lty=2,col="steelblue",lwd=1)
      title(sub=paste("Scenario: delta T=",round(dT,1),"C, ", 
            " delta P=",round(dP,1),"mm/day",sep=""))
      text(70,0.08,model.chg,cex=0.8,col="blue")
    }
  }  

  cdf.obs <- cumsum(pdf)/sum(pdf)
  cdf.chg <- cumsum(pdf.chg)/sum(pdf.chg)
 
  results <- list(fx.obs=pdf,fx.chg=pdf.chg,x=exp.x,Fx.obs=cdf.obs,Fx.chg=cdf.chg,
                  model=model,location=obs$location,lon=obs$lon,lat=obs$lat,alt=obs$alt,
                  minAmountPrecip=exp.par$minAmountPrecip,dT=dT,dP=dP,
                  slope.dep=slope.dep, slope.x=slope.chg,slope.coef=summary(smod)$coefficients)
  invisible(results)
}


CDFtransfer <-  function(Y,CDF.2,CDF.1=NULL,method="empiricalRanking",
                         plot=FALSE,silent=FALSE,smooth=TRUE) {

  if (class(Y)[1]=="station") {
    obs <- Y
    print("Extracting precip from station object")
    Y <- obs$precip
  }
  
  if (is.null(CDF.1)) {
    if (!silent) print("Using the emprical distribution function")
    CDF.1 <- eval(parse(text=paste(method,"(Y)",sep="")))
  }

  minmax <- range(c(0,CDF.1$x,CDF.2$x),na.rm=TRUE)
  F1 <- spline(CDF.1$x,CDF.1$P,n=300,xmin=minmax[1],xmax=minmax[2])
  F2 <- spline(CDF.2$x,CDF.2$P,n=300,xmin=minmax[1],xmax=minmax[2])

  x1 <- seq(minmax[1],minmax[2],length=100)
  x2 <- rep(NA,length(x1)); prob <- x2
  
  for (i in 1:length(x1)) {
    i1 <- (F1$x <= x1[i])
    prob[i] <- max(F1$y[i1],na.rm=TRUE)
    i2 <- (F2$y <= prob[i])
    x2.mn <- max(F2$x[i2],na.rm=TRUE)

    i1 <- (F1$x >= x1[i])
    prob[i] <- min(F1$y[i1],na.rm=TRUE)
    i2 <- (F2$y >= prob[i])
    x2mx <- min(F2$x[i2],na.rm=TRUE)

    if ( (sum(i1)>0) & (sum(i2)>0) ) x2[i] <- mean(c(x2.mn,x2mx)) else {
                                     x2[i] <- minmax[2]
      print(paste("CDFtransfer",i,x1[i],x2.mn,x2mx,x2[i],sum(i1),sum(i2)))
    }
  }

  s <- 2*sd(x1,na.rm=TRUE); m <- mean(x1,na.rm=TRUE)
  x1x2 <- data.frame(y=(x2-m)/s,x=(x1-m)/s)
  x.smooth <- lm(y ~ x+I(x^2)+ I(x^3)+I(x^4)+I(x^5)+I(x^6)+I(x^7)+I(x^8)+I(x^9),data=x1x2)
  
  if (plot) {
    x11()
    plot(x1,x2,main="Local quantile transfer function",type="n")
    grid()
    lines(minmax,minmax,col="grey70")
    points(x1,x2,pch=20,col="grey30",cex=0.7)
    lines(x1,predict(x.smooth)*s + m,col="red")
  }  
  if (smooth) x2 <- predict(x.smooth)*s+m
  
  Y.new <- rep(NA,length(Y))
  for (i in 1:length(Y)) Y.new[i] <- min(x2[(x1 >= Y[i])],na.rm=TRUE)
  
  if (exists("obs",envir=environment(CDFtransfer))) {
    obs$precip <- Y.new
    Y.new <- obs
  }
  invisible(Y.new)
}

# Empirical ranking method:
#
# A formula for estimating the cumulative probability P corresponding to rank m
# Original references:
#  Jenkinson, A.F., 1977, U.K. Met.Office Synoptic Clim. Branch Memo 58
#  Beard, L.R., 1943, Trans. Amer. Meteor. Soc. Civ. Eng., 108, 1110-1160
#  Chegodaev, N.N., 1953 (in Russian) State Rail Transport Publishing House.
# Reference:
# Folland, C. and Anderson, C. (2002), J. Clim. 15, 2954-2960, equation (1)
#
empiricalRanking <- function(x) {
  N <- length(x)
  m <- rank(x)
  P <- (m-0.31)/(N+0.38)
  sort <- order(x)
  P <- P[sort]
  y <- as.numeric(x[sort])
  results <- data.frame(x=y, P=P)
  results
}
