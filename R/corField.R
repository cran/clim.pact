corField <- function(x,y,lsig.mask=TRUE,sig.lev=0.05,mon=NULL,param="t2m",
                     lty=1,col="black",lwd=1,main=NULL,z.levs=NULL,my.col=NULL,plot=TRUE) {
#  library(ctest)
  library(akima)
  
if ((class(x)[1]!="field") & (class(x)[1]!="monthly.field.object") &
    (class(x)[1]!="daily.field.object")){
  stop("x must be a 'field' object.")
}

if ((class(y)[1]!="station") & (class(y)[1]!="field") &
    (class(y)[1]!="monthly.field.object")) {
  stop(paste("y must be a 'monthly.station.record' or a 'field'",
             "object - Use  station.obj()"))
}

if (is.null(mon)) mon <- 1:12
  
 cmon<-c('Jan','Feb','Mar','Apr','May','Jun',
         'Jul','Aug','Sep','Oct','Nov','Dec') 
  descr <- 'Correlation:'
  date <- ""

  if (class(x)[2]=="monthly.field.object") {
    d <- dim(x$dat); dim(x$dat) <- c(d[1],d[2]*d[3])                                 # REB 12.01.2006
    nyears <- max(x$yy)-min(x$yy)+1                                                  # REB 12.01.2006
    dat <- rep(NA,nyears*d[2]*d[3]); dim(dat) <- c(nyears,d[2]*d[3])                 # REB 12.01.2006
    yy <- rep(NA,nyears); mm <- yy; dd <- yy; tim <- yy; id.t <- yy                  # REB 12.01.2006
    if (!is.null(mon)) {
      iyear <- 0                                                                     # REB 12.01.2006
      for (year in min(x$yy):max(x$yy)) {                                            # REB 12.01.2006
        iyear <-iyear + 1                                                            # REB 12.01.2006
# REB before 12.01.2006:         im <- is.element(x$mm,mon) 
        im <- is.element(x$mm,mon) & is.element(x$yy,year)                           # REB 12.01.2006
        #print(c(iyear,year,sum(im)))
# REB before 12.01.2006:   x$dat <- x$dat[im,,]
        if (sum(im)>1) dat[iyear,] <- colMeans(x$dat[im,]) else                      # REB 12.01.2006
          if (sum(im)==1) dat[iyear,] <- x$dat[im,]                                  # REB 12.01.2006
        if (sum(im) > 0) {
          yy[iyear] <- x$yy[im]                                                      # REB 12.01.2006
          mm[iyear] <- mon[1]; dd[iyear] <- 15                                       # REB 12.01.2006
          tim <- x$tim[im]                                                           # REB 12.01.2006
          id.t <- x$id.t[im]                                                         # REB 12.01.2006
        }
      }                                                                              # REB 12.01.2006
      #print(c(dim(x$dat),NA,d,NA,nyears))
      dim(dat) <- c(nyears,d[2],d[3])                                        #bug-fix  REB 06.08.2007
# REB 12.01.2006
      if (length(mon)==1) date <- cmon[mon] else
                          date <- paste(cmon[mon[1]],"-",cmon[mon[length(mon)]])
      date <- paste(date,": ",min(x$yy)," - ",max(x$yy),sep="")
      } else {
        nyears <- d[1]                                                           # REB 31.08.2007
        dim(x$dat) <- c(d[1],d[2],d[3])                                          # REB 31.08.2007
        dat <- x$dat; yy <- x$yy; mm <- x$mm; tim <- x$tim; id.t <- x$id.t       # REB 31.08.2007    
      }                                                                          # endif (!is.null(mon))
  } else {                                                                           # REB 05.07.2007
      if (!is.null(mon)) ii <- is.element(x$mm,mon) else
                         ii <- is.finite(x$mm)
      if (length(mon)==1) date <- cmon[mon] else
                          date <- paste(cmon[mon[1]],"-",cmon[mon[length(mon)]])
      yy <- x$yy[ii]; mm <- x$mm[ii]; dd <- x$dd[ii]; id.t <- x$id.t[ii]
      dat <- x$dat[ii,,]; tim <- x$tim[ii]
  }
  ok <- is.finite(yy)
  x$dat <- dat[ok,,]; x$yy <- yy[ok]; x$mm <- mm[ok]; x$dd <- dd[ok];
  x$tim <- tim[ok]; x$id.t <- id.t[ok]                                              # REB 12.01.2006
  nyears <- sum(ok)
  rm(dat,yy,mm,dd,tim,id.t)                                                         # REB 12.01.2006

# REB before 12.01.2006:  if (class(y)[1]=="station") y.ts <- as.vector(t(y$val)) else
  if (class(y)[1]=="station") {
    if (class(y)[2]=="monthly.station.record") {
#    if (length(mon) > 1) y.ts <- rowMeans(y$val[,mon]) else                           # REB 12.01.2006
                         y.ts <- y$val[,mon]                                          # REB 12.01.2006
    param <- y$obs.name                                                               # REB 05.07.2007
  } else if (class(y)[2]=="daily.station.record") {                                   # REB 05.07.2007
    y.ts <- eval(parse(text=paste("y$",param,sep="")))                                # REB 05.07.2007
  }
  #print(paste("length(y.ts)=",length(y.ts),"dim(y$val[,mon])=",dim(y$val[,mon])))
  } else {
    l.diffgrid <- TRUE
    if ( (length(x$lon)==length(y$lon)) & (length(x$lat)==length(y$lat)) ) {
      if ( (sum(x$lon==y$lon)==length(x$lon)) &
           (sum(x$lat==y$lat)==length(x$lat)) ) l.diffgrid <- FALSE
    }
    y.ts <- rep(NA,length(y$tim)*length(x$lat)*length(x$lon))
    dim(y.ts) <- c(length(y$tim),length(x$lat),length(x$lon))
    if (l.diffgrid) {
      print("Interpolate 2nd field - please be patient :-)")
      lat.x<-rep(x$lat,length(x$lon))
      lon.x<-sort(rep(x$lon,length(x$lat)))
      for (i in 1:length(y$tim)) {
         map <- as.matrix(y$dat[i,,])
         y.ts[i,,] <- interp(lat.x,lon.x,map,y$lat,y$lon)$z
      }
    } else {
      for (i in 1:length(y$tim)) {
         y.ts[i,,] <- as.matrix(y$dat[i,,])
      }
    }
  }
    
  if (class(y)[1]=="station") {
    good <- is.finite(y.ts)
# REB before 12.01.2006:       yy <- sort(rep(y$yy,12))
# REB before 12.01.2006:       mm <- rep(1:12,length(y$yy))
# REB before 12.01.2006:       dd <- rep(15,length(yy))
    yy <- y$yy
    if (class(y)[2]=="monthly.station.record") {                                     # REB 05.07.2007
      mm <- rep(mon,length(yy))
      if (length(mon)>1) yy <- sort(rep(yy,length(mon)))
      dd <- rep(15,length(yy))
#      print(summary(mm)); print(class(mm)); print(length(yy))
    } else {mm <- y$mm; dd <- y$dd}                                                # REB 05.07.2007
    #print(sum(good))
    y.ts <- y.ts[good]; yy <- yy[good]; mm <- mm[good]; dd <- dd[good]
  } else {
    yy <- y$yy; mm <- y$mm;  dd <- y$dd
  }
  
    i1<-is.element(yy*10000+mm*100+dd,
                   x$yy*10000+x$mm*100+x$dd)
    i2<-is.element(x$yy*10000+x$mm*100+x$dd,
                   yy*10000+mm*100+dd)

  if ( (sum(i1)==0) | (sum(i2)==0) ) {
    print("Years:")
    print(table(yy))
    print(table(x$yy))
    print("Months:")
    print(table(mm))
    print(table(x$mm))
    print("(days):")
    print(table(dd))
    print(table(x$dd))
    print(c(sum(i1),sum(i2)))
    stop('CorField: no matching dates.')
  }
  #print(rbind(yy[i1]*100+mm[i1],x$yy[i2]*100+x$mm[i2]))
  ni <- length(x$lon)
  nj <- length(x$lat)
  map <- matrix(rep(NA,ni*nj),nj,ni)
  p.val <- matrix(rep(NA,ni*nj),nj,ni)

#  print("yy:"); print(range(yy)); print(range(x$yy)); print("mm:"); print(range(mm)); print(range(x$mm))
#  print("dd:"); print(range(dd)); print(range(x$dd)); print("i1/i2:"); print(c(sum(i1),sum(i2)))
#  print("dim:"); print(dim(x$dat)); print(dim(x$dat[i2,,]));
#  print(dim(y.ts)); print(dim(y.ts[i1,,])); print(class(y)[1])
#  print(table(yy))
#  print(dim(x$dat[i2,,]));  print(length(y.ts)); print(length(y.ts[i1])); print(class(y)[1])

    for (j in 1:nj) {
      for (i in 1:ni) {

        if ( (class(y)[1]=="station") & (sum(is.finite(x$dat[i2,j,i]))> 10) ) {
                                    good <- is.finite(y.ts[i1]) & is.finite(x$dat[i2,j,i])
#                                    ii1 <- i1 & good
#                                    ii2 <- i2 & good
#                                    ii1 <- is.finite(y.ts[i1])
#                                    ii2 <- is.finite(x$dat[i2,j,i])
#        if (length(i1) != length(ii1)) print(c(length(i1),length(ii1),length(good),
#                                               length(y.ts[i1]),length(x$dat[i2,j,i])))
#        if (length(i2) != length(ii2)) print(c(length(i2),length(ii2),length(good),
#                                               length(y.ts[i1]),length(x$dat[i2,j,i])))
#        print(rbind(yy[ii1]*100+mm[ii1],
#                    x$yy[ii2]*100+x$mm[ii2]))
                                    r.test <- cor.test(x$dat[i2,j,i][good],y.ts[i1][good])
                                  } else if (class(y)[1]=="field") {
                                    good <- is.finite(y.ts[i1,j,i]) & is.finite(x$dat[i2,j,i])
#                                    ii1 <- is.finite(y.ts[i1,j,i])
#                                    ii2 <- is.finite(x$dat[i2,j,i]
                                    r.test <- cor.test(x$dat[i2,j,i][good],y.ts[i1,j,i][good])
                                  } else {
                                    #print(c(i,j,sum(is.finite(x$dat[i2,j,i]))))
                                    r.test <- list(estimate=NA,p.value=NA) }
        map[j,i] <- r.test$estimate
        p.val[j,i] <- r.test$p.value
      }
      #good <- is.finite(y.ts[i1]) & is.finite(x$dat[i2,j,i]); ii1 <- i1 & good; ii2 <- i2 & good
      #print(c(sum(i1),sum(i2),sum(ii1),sum(ii2),sum(is.finite(x$dat[,j,i])),sum(is.finite(x$dat[i2,j,i])),
      #        sum(is.finite(y.ts[i1])),map[j,i],p.val[j,i]))
    }
  if (is.null(z.levs)) {
    z.levs <- seq(-1,1,length=41)
  }
  if (is.null(my.col)) {
    my.col <- rgb(c(seq(0,1,length=20),rep(1,21)),
                  c(abs(sin((0:40)*pi/40))),
                  c(c(rep(1,21),seq(1,0,length=20))))
  }

#  print(dim(map))
#  print(c(length(x$lon),length(x$lat)))

  if (is.null(attributes(x$dat)$"long_name")) attr(x$dat,"long_name") <- x$v.name
  if (is.null(attributes(y$dat)$"long_name")) attr(y$dat,"long_name") <- y$v.name
  
  if (is.null(main)) {
    if (class(y)[1]=="station") {
      main <- paste(descr,x$v.name,"&",param,"at",y$location)
      if (!is.null(attributes(y)$lagStation)) date <- paste(date," (lag=",attributes(y)$lagStation,")",sep="")
    } else main <- paste(descr,x$v.name,"&",y$v.name,"(Pointwise)")
    descr <- main
  }

  if (lsig.mask) {
     if (sum(p.val < 0.05,na.rm=TRUE) > 0)  map[p.val > 0.05] <- NA else {
      my.col <- rgb(rep(1,21),rep(1,21),rep(1,21)); col="grey80"}
  }
  

  if (plot) {
  newFig()
  filled.contour(x$lon,x$lat,t(map),
                 col = my.col,levels=z.levs,
                 main=main,
                 sub=date,xlab="Longitude",ylab="Latitude")

# From filled.contour in base
  mar.orig <- (par.orig <- par(c("mar","las","mfrow")))$mar
  on.exit(par(par.orig))

  w <- (3 + mar.orig[2]) * par('csi') * 2.54
  layout(matrix(c(2, 1), nc=2), widths=c(1, lcm(w)))
    
  par(las = 1)
  mar <- mar.orig
  mar[4] <- 1
  par(mar=mar)
  contour(x$lon,x$lat,t(map),add=TRUE,col=col,lwd=lwd,lty=lty)
  addland()
  if (class(y)[1]=="station") {
    points(y$lon,y$lat,pch=20,col="white",cex=1.2)
    points(y$lon,y$lat,pch=20,col="black",cex=0.9)
  }
  }
  results <- list(map=t(map),lon=x$lon,lat=x$lat,tim=x$tim,
                  date=date,description=descr)
  class(results) <- "map"
  attr(results,"long_name") <- attr(x$dat,"long_name")
  attr(results,"descr") <- "Correlation map"
  invisible(results)
}
