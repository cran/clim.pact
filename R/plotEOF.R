# Plots EOF products
# Monthly mean values.
#
# Reference: R.E. Benestad et al. (2002),
#            Empirically downscaled temperature scenarios for Svalbard,
#            submitted to Atm. Sci. Lett.
#
#            R.E. Benestad (2001),
#            A comparison between two empirical downscaling strategies,
#            Int. J. Climatology, 1645-1668, vol. 21, DOI 10.1002/joc.703
#
# R.E. Benestad, met.no, Oslo, Norway 10.05.2002
# rasmus.benestad@met.no
#

#------------------------------------------------------------------------

plotEOF<-function(x,i.eof=1,nlevs=5,
                   col=c("red","blue","darkgreen","steelblue"),
                   main=NULL,sub=NULL) {

if (class(x)[1]!= "eof") stop ("The argument must be an 'eof' object") 
ok.eps <- (lower.case(options()$device)=="x11") | (lower.case(options()$device)=="windows")
attach(x)

dims <- dim(x$EOF) 
if (length(dims)==3) dim(EOF) <- c(dims[1],dims[2]*dims[3])

title.1 <- paste("EOF pattern #",i.eof,"(",class(x)[2],")",sep="")
title.2 <- "The fraction of variance accounted by the EOFs"
title.3 <- paste("Principal component (",class(x)[2],")",sep="")

vnames <- x$v.name[1]
if (!is.null(main)) {
   title.1  <-  main; title.2 <- main; title.3 <- main}       
if (is.null(sub)) sub <- paste(vnames," (",c.mon,")")

for (i in 2:length(vnames)) vnames <- paste(vnames,"+",x$v.name[i])
i.last <- 0
id <- row.names(table(id.x))
ordr <- rep(NA,length(id))
for (i in 1:length(id)) {
  ordr[i] <- min((1:length(id.x))[is.element(id.x,id[i])])
}
#print(ordr)
id<-id[order(ordr)]
#print(id)

#par(ask=TRUE)
newFig()
plot(c(floor(min(x$lon)),ceiling(max(x$lon))),
     c(floor(min(x$lat)),ceiling(max(x$lat))),
     type="n",main=title.1,sub=sub,
     xlab="Longitude",ylab="Latitude")
if (range(x$lon)[2]-range(x$lon)[1] > 360) {
  xy.cont <- COn0E65N(lon.cont, lat.cont)
  addland(lon=xy.cont$x,lat=xy.cont$y)
} else addland()
grid()
col.tab <- col[1:length(id)]
neofs <- length(x$var)
i.last <- 0
print(paste("plotEOF: n.fld=",n.fld))
for (i in 1:n.fld) {
  i.fld <- seq(i.last+1,i.last+x$size[2,i]*x$size[3,i],by=1)
  i.last <- max(i.fld)
  #print(c(i,NA,dim(x$EOF),NA,size[,i],NA,range(i.fld)))
  EOF.1 <- x$EOF[,i.fld]
  dim(EOF.1)<-c(dim(x$EOF)[1],size[2,i],size[3,i])
  eof.patt<-t(EOF.1[i.eof,,])
  i.lon <- x$id.lon == id[i]
  i.lat <- x$id.lat == id[i]
  lon.x <- x$lon[i.lon]
  lat.x <- x$lat[i.lat]
  #print(c(size[,i],NA,length(lon.x),length(lat.x),NA,dim(eof.patt),id[i]))
  contour(lon.x,lat.x,eof.patt,
          nlevels=nlevs,add=TRUE,lwd=2,col=col.tab[i])
}

if (n.fld>1) legend(min(x$lon),max(x$lat),id,
             col=c(col.tab),lty=1,
             lwd=2,merge=TRUE, bg='gray95')

if (ok.eps) dev.copy2eps(file=paste("plotEOF_1.eps",sep=""))

newFig()
plot(100*(W+dW)^2/tot.var,main=title.2,type="n",
     ylab="Variance (%)",xlab="EOF order",sub=sub)
lines(var.eof,lty=3)
for (i in 1:length(var.eof)) {
  lines(rep(i,2),100*c((W[i]+dW[i])^2/tot.var,(W[i]-dW[i])^2/tot.var),
        lty=2,col="darkgrey")
  lines(c(i-0.25,i+0.25),100*rep((W[i]+dW[i])^2/tot.var,2),
        lwd=2,col="darkgrey")
  lines(c(i-0.25,i+0.25),100*rep((W[i]-dW[i])^2/tot.var,2),
        lwd=2,col="darkgrey")
}
points(var.eof)
points(var.eof,pch=20,cex=0.8,col="darkgrey")
grid()
if (ok.eps) dev.copy2eps(file=paste("plotEOF_2.eps",sep=""))

newFig()
yymm<-x$yy + (x$mm-0.5)/12 + (x$dd-0.5)/365.25
#print(c(length(yy),length(mm),length(dd),length(yymm),length(PC[,i.eof])))
plot(yymm,PC[,i.eof],pch=20,cex=0.7,
     main=title.3,,col="grey70",sub=sub)

lines(yymm[id.t==id.t[1]],PC[id.t==id.t[1],i.eof],col="red",lty=2,lwd=2)
if (sum(id.t!=id.t[1])>0) lines(yymm[id.t!=id.t[1]],
          PC[id.t!=id.t[1],i.eof],col="blue",lty=2,lwd=2)
grid()
if (ok.eps) dev.copy2eps(file=paste("plotEOF_3.eps",sep=""))

detach(x)
}

