# Plots the spatial EOF patterns 
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


map.eof <- function(x,i.eof=1,nlevs=5,add=FALSE,
            col=c("red","blue","darkgreen","steelblue"),lwd=2,lty=1) {

if (class(x)!= "eof") stop ("The argument must be an 'eof' object") 

attach(x)

dims <- dim(EOF) 
if (length(dims)==3) dim(EOF) <- c(dims[1],dims[2]*dims[3])

title.1 <- paste("EOF pattern #",i.eof,"(",class(x)[2],")",sep="")
i.last <- 0
id <- row.names(table(id.x))
if (!add) {
  plot(c(floor(min(lon)),ceiling(max(lon))),
     c(floor(min(lat)),ceiling(max(lat))),
     type="n",main=title.1,
     sub=paste(x$f.name," (",c.mon,")"),
     xlab="Longitude",ylab="Latitude")
}

if (range(lon)[2]-range(lon)[1] > 360) {
  xy.cont <- COn0E65N(lon.cont, lat.cont)
  addland(lon=xy.cont$x,lat=xy.cont$y)
} else addland()
grid()

col.tab <- col[1:length(id)]
neofs <- length(x$var)
i.last <- 0
for (i in 1:n.fld) {
  i.fld <- seq(i.last+1,i.last+size[2,i]*size[3,i],by=1)
  i.last <- max(i.fld)
  EOF.1 <- EOF[,i.fld]
  dim(EOF.1)<-c(neofs,size[2,i],size[3,i])
  eof.patt<-t(EOF.1[i.eof,,])
  i.lon <- id.lon == id[i]
  i.lat <- id.lat == id[i]
  lon.x <- lon[i.lon]
  lat.x <- lat[i.lat]
  contour(lon.x,lat.x,eof.patt,
          nlevels=nlevs,add=TRUE,lwd=2,col=col.tab[i])
}
if ((n.fld>1) & (!add)) legend(min(lon),max(lat),id,
             col=c(col.tab),lty=1,
             lwd=2,merge=TRUE, bg='gray95')

  cmon<-c('Jan','Feb','Mar','Apr','May','Jun',
          'Jul','Aug','Sep','Oct','Nov','Dec')

  results <- list(map=eof.patt,lon=x$lon,lat=x$lat,tim=x$tim,
                  date=c.mon,description=id[i])

  class(results) <- c("map","eof")
}
