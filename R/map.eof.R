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

map.eof <- function(x,i.eof=1,nlevs=9,add=FALSE,
            col=c("red","blue","darkgreen","steelblue"),lwd=2,lty=1) {

if (class(x)[1]!= "eof") stop ("The argument must be an 'eof' object") 

attach(x)

dims <- dim(x$EOF) 
if (length(dims)==3) dim(x$EOF) <- c(dims[1],dims[2]*dims[3])

title.1 <- paste("EOF pattern #",i.eof,"(",class(x)[2],")",sep="")
i.last <- 0
id <- row.names(table(x$id.x))
if (!add) {
  plot(c(floor(min(lon)),ceiling(max(x$lon))),
     c(floor(min(lat)),ceiling(max(x$lat))),
     type="n",main=title.1,
     sub=paste(x$f.name," (",c.mon,")"),
     xlab="Longitude",ylab="Latitude")
}

#if (range(x$lon)[2]-range(x$lon)[1] > 360) {
#  xy.cont <- COn0E65N(lon.cont, lat.cont)
#  addland(lon=xy.cont$x,lat=xy.cont$y)
#} else
addland()
grid()

col.tab <- col[1:length(id)]
neofs <- length(x$var)
i.last <- 0
for (i in 1:x$n.fld) {
  #print(c(i,NA,x$size[,i],NA,i.last,NA,i.last+1,i.last+x$size[2,i]*x$size[3,i]))
  i.fld <- seq(i.last+1,i.last+x$size[2,i]*x$size[3,i],by=1)
  i.last <- max(i.fld)
  EOF.1 <- x$EOF[,i.fld]
  dim(EOF.1)<-c(neofs,x$size[2,i],x$size[3,i])
  eof.patt<-t(EOF.1[i.eof,,])
  i.lon <- x$id.lon == id[i]
  i.lat <- x$id.lat == id[i]
  lon.x <- x$lon[i.lon]
  lat.x <- x$lat[i.lat]
  #print(summary(as.vector(eof.patt)))
  #print(lon.x)
  #print(lat.x)
  contour(lon.x,lat.x,eof.patt,
          nlevels=nlevs,add=TRUE,lwd=2,col=col.tab[i])
}
if ((x$n.fld>1) & (!add)) legend(min(x$lon),max(x$lat),id,
             col=c(col.tab),lty=1,
             lwd=2,merge=TRUE, bg='gray95')

  cmon<-c('Jan','Feb','Mar','Apr','May','Jun',
          'Jul','Aug','Sep','Oct','Nov','Dec')

  results <- list(map=eof.patt,lon=x$lon,lat=x$lat,tim=x$tim,
                  date=c.mon,description=id[i])

  class(results) <- c("map","eof")
}


