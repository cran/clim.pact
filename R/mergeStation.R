# Merges two station series of the same variable but from different sources to
# produce a long, updated series.
# R.E. Benestad

mergeStation <- function(x.1,x.2,plot=FALSE,print=TRUE,rescale=TRUE) {

if ( (class(x.1)[2]!="monthly.station.record") |
     (class(x.2)[2]!="monthly.station.record")) {
  stop(paste("The predictand must be a 'monthly.station.record'",
             "object - Use  station.obj()"))
}
#print("Time intervals:")
#print(range(x.1$yy))
#print(range(x.2$yy))

if (min(x.1$yy) > min(x.2$yy)) {
  XX <- x.1
  x.1 <-  x.2
  x.2 <- XX
  rm(XX)
}
ny.1 <- length(x.1$yy)
ny.2 <- length(x.2$yy)
y.1 <- t(x.1$val)
y.2 <- t(x.2$val)
dim(y.1) <- c(12*ny.1); y.1 <- as.vector(y.1)
dim(y.2) <- c(12*ny.2); y.2 <- as.vector(y.2)
y.1[y.1 <= -999] <- NA
y.2[y.2 <= -999] <- NA
yymm.1 <- sort(rep(x.1$yy,12)) + (rep(1:12,ny.1)-0.5)/12
yymm.2 <- sort(rep(x.2$yy,12)) + (rep(1:12,ny.2)-0.5)/12
i.1 <- is.element(yymm.1,yymm.2)
i.2 <- is.element(yymm.2,yymm.1)
#print(yymm.1); print(i.1)
if (sum(i.1)) {
  if (print) print(range(yymm.1[i.1]))
  ovrlp <- data.frame(y=y.1[i.1],x=y.2[i.2])
  new.dat <- data.frame(x=y.2[!i.2])

  Y.1 <- y.1[i.1]
  Y.2 <- y.2[i.2]
  ii <- is.finite(Y.1) & is.finite(Y.2)

  if (print) print(paste("RMSE: ",round(sqrt(sum( (Y.1[ii]-Y.2[ii])^2 ))/sum(i.1),2)))

  agree <- lm(y ~ 1 + x, data=ovrlp)
  if (print) print(summary(agree))
  coefs <- agree$coefficients
  if (!rescale) coefs[2] <- 1

#print("New series")
  y <- c(y.1,coefs[1] + coefs[2]* y.2[!i.2])
#print("Years")
  yy <- c(x.1$yy,as.numeric(row.names(table(floor(yymm.2[!i.2])))))
  ny <- length(yy)
  yymm <- sort(rep(yy,12)) + (rep(1:12,ny) - 0.5)/12
} else {
   y <- c(y.1,y.2)
   yy <- c(x.1$yy,x.2$yy)
   #print(c(length(yy),length(y),NA,length(yy)*12))
   ny <- length(yy)
   yymm <- sort(rep(yy,12)) + (rep(1:12,ny) - 0.5)/12
}
#print("Plot?")
if (plot) {
  plot(yymm.1,y.1,type="s",lwd=3,col="darkblue",xlim=range(yymm))
  lines(yymm.2,y.2,type="s",col="steelblue",lty=3,lwd=2)
  lines(yymm,y,type="s",col="wheat",lwd=1)
  grid()
}

#print("Change dimensions")
dim(y) <- c(12,ny)

x<- x.1
x$val <- t(y)
x$yy <- yy
x
}

