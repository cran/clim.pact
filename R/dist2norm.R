# Transform to normal distribution 

dist2norm <- function(x,plot=FALSE,exclude=NULL,
                      sd=1,mean=0,force.zero=TRUE) {
  xT <- rep(NA,length(x))
  good <- is.finite(x)
  if (!is.null(exclude)) good <- good & (x!=exclude)
  breaks <- seq(min(x[good])-0.5*IQR(x[good]),
                max(x[good])+0.5*IQR(x[good]),length=15)
  h <- hist(x[good],breaks=breaks,plot=plot)
  q <- h$mids
  q.i <- seq(-3*sd+mean,3*sd+mean,length=100)
  edf <- cumsum(h$density)/sum(h$density)
  if (force.zero) {
    edf <- (edf-edf[1])
    edf <- edf/max(edf)
  }
  edf.i <- spline(q,edf,n=100)
  cdf.i <- pnorm(q.i,mean=mean,sd=sd)
  if (plot) {
    x11()
    plot(c(min(c(q.i,edf.i$x)),max(c(q.i,edf.i$x))),c(0,1),type="n",
         main="dist2normal",xlab="Value",ylab="Probability")
    grid()
    lines(edf.i$x,edf.i$y,type="l",lwd=2)
    lines(q.i,cdf.i,col="red")
  }

    for (i in 1:length(x)) {
    if (is.finite(x[i])) {
      y.i <- approx(edf.i$x,edf.i$y,x[i])$y
      xT[i] <- approx(cdf.i,q.i,y.i)$y
#      print(c(i,round(x[i],2),round(xT[i],2),round(y.i,2)))
      if ((mod(i,round(length(x)/10))==0) & (plot)){
        lines(rep(x[i],2),c(0,y.i),col="blue",lty=2)
        lines(c(x[i],xT[i]),rep(y.i,2),col="blue",lty=2)
        arrows(xT[i],y.i,xT[i],0,col="blue",lty=2,length=0.05)
      }
    }
  }

  dist2norm <- list(xT=xT,x=x,edf.i=edf.i,cdf.i=cdf.i,q.i=q.i)
  class(dist2norm) <- "dist2norm"
  invisible(dist2norm)
}


norm2dist <- function(x,plot=FALSE,exclude=NULL,
                      sd=1,mean=0,force.zero=TRUE) {
  if (class(x)!="dist2norm") {
    print("Needs a 'dist2norm'-type object")
    return()
  }
  if (plot) {
    x11()
    plot(c(min(c(x$q.i,x$edf.i$x)),max(c(x$q.i,x$edf.i$x))),c(0,1),type="n",
         main="dist2normal",xlab="Value",ylab="Probability")
    grid()
    lines(x$edf.i$x,x$edf.i$y,type="l",lwd=2)
    lines(x$q.i,x$cdf.i,col="red")
  }
  for (i in 1:length(x$xT)) {
    if (is.finite(x$xT[i]))
      print(approx(x=x$q.i,y=x$cdf.i,xout=x$xT[i]))$y
      Y.i <- approx(x=x$q.i,y=x$cdf.i,xout=x$xT[i])$y
      if (is.finite(Y.i)) {
        x$x[i] <- approx(x=x$edf.i$x,y=x$cdf.i,xout=Y.i)$y
        if ((mod(i,round(length(x$x)/10))==0) & (plot)){
          lines(rep(x$xT[i],2),c(0,Y.i),col="blue",lty=2)
          lines(c(x$x[i],x$xT[i]),rep(Y.i,2),col="blue",lty=2)
          arrows(x$x[i],Y.i,x$x[i],0,col="blue",lty=2,length=0.05)
        }
      }
  }
  norm2dist <- x$x
  invisible(norm2dist)
}
