CCA <- function(x1,x2,SVD=TRUE,plot=TRUE,main="CCA",sub="",test=FALSE,i.eofs=1:8,LINPACK=TRUE) {

  i1 <- is.element(x1$yy*10000 + x1$mm*100 + x1$dd, x2$yy*10000 + x2$mm*100 + x2$dd)
  i2 <- is.element(x2$yy*10000 + x2$mm*100 + x2$dd, x1$yy*10000 + x1$mm*100 + x1$dd)

  if ( (class(x1)[1]=="field") & (class(x2)[1]=="field") |
       (class(x1)[1]=="monthly.field.object") & (class(x2)[1]=="monthly.field.object") ) {
    print("classical CCA")
    X1 <- x1$dat[i1,,]; X2 <- x2$dat[i2,,]
    d.1 <-dim(X1); d.2 <- dim(X2)
    dim(X1) <- c(d.1[1],d.1[2]*d.1[3])
    dim(X2) <- c(d.2[1],d.2[2]*d.2[3])
    X1.m <- rowMeans(X1); X2.m <- rowMeans(X2)
    X1 <- X1 - X1.m; X2 <- X2 - X2.m
    C.12 <- cov(X1,X2)
    C.11 <- cov(X1,X1)
    C.22 <- cov(X2,X2)

    if (SVD) {
# CCA After Bretherton et al. (1992), J. Clim. Vol 5, p. 541:
      d1 <- dim(C.11); d2 <- dim(C.22)
      C.11 <- round(as.complex(C.11),4); dim(C.11) <- d1
      C.22 <- round(as.complex(C.22),4); dim(C.22) <- d2
      C <- Re(sqrt(solve(C.11)) %*% C.12 %*%  sqrt(solve(C.22)))
      if (LINPACK) M <-svd(C) else 
                   M <-La.svd(C)
      l.k <- M$u; r.k <- M$v; R <- M$d
      a.m <- l.k %*% C.11 
      b.m <- r.k %*% C.22 
      u.k <- X1 %*% t(a.m)
      v.k <- X2 %*% t(b.m)
    } else {

# After Wilks, 1995, p. 401
      sub <- paste(sub,"(BP CCA - after Wilks (1995))")
      M.1 <- solve(C.11) %*% C.12 %*% solve(C.22) %*% C.12
      M.2 <- solve(C.22) %*% C.12 %*% solve(C.11) %*% C.12
      a.m <- eigen(M.1)
      b.m <- eigen(M.2)
print(dim(a.m)); print(dim(X1))
      u.k <- X1 %*% a.m
      v.k <- X2 %*% b.m
      R <- sqrt(Re(x.m$values))
   }
    dim(a.m) <- c(d1[1],d.1[2],d.1[3]); dim(b.m) <- c(d2[1],d.2[2],d.2[3])


  } else if ( (class(x1)[1]=="eof") & (class(x2)[1]=="eof") ) {

    n.eof1 <- length(x1$PC[1,]);   n.eof2 <- length(x2$PC[1,])
    i.eofs <- i.eofs[(i.eofs <= n.eof1) & (i.eofs <= n.eof2)]

    X1 <- x1$PC[i1,i.eofs]; X2 <- x2$PC[i2,i.eofs]
    print("Barnett-Preisendorfer CCA")
    C.12 <- cov(X1,X2)
    C.11 <- cov(X1,X1)
    C.22 <- cov(X2,X2)

    d.1 <-c(length(x1$W[i.eofs]),length(x1$lat),length(x1$lon))
    d.2 <-c(length(x2$W[i.eofs]),length(x2$lat),length(x2$lon))
    if (SVD) {
      sub <- paste(sub,"(BP CCA - after Bretherton et al. (1992))")
      print("This method gives somewhat strange results - do not trust it!")
# CCA After Bretherton et al. (1992), J. Clim. Vol 5, p. 541:
      d1 <- dim(C.11); d2 <- dim(C.22)
# C.11 & C.22 should be nearly diagonal.
      C <- diag(1/sqrt(diag(C.11))) %*% C.12 %*%  diag(1/sqrt(diag(C.22)))
     if (LINPACK) M <-svd(C) else 
                  M <-La.svd(C)
      l.k <- M$u; r.k <- M$v; R <- M$d
      a.m <- t( t(x1$EOF[i.eofs,]) %*% diag(x1$W[i.eofs]) %*% l.k )
      b.m <- t( t(x2$EOF[i.eofs,]) %*% diag(x2$W[i.eofs]) %*% r.k )
      u.k <- x1$PC[,i.eofs] %*% l.k
      v.k <- x2$PC[,i.eofs] %*% r.k

    } else {

# After Wilks, 1995, p. 401
      sub <- paste(sub,"(BP CCA - after Wilks (1995))")
      M.1 <- solve(C.11) %*% C.12 %*% solve(C.22) %*% C.12
      M.2 <- solve(C.22) %*% C.12 %*% solve(C.11) %*% C.12
      x.m <- eigen(M.1)
      y.m <- eigen(M.2)

      #print(dim( t(x1$EOF[i.eofs,]) )); print(dim( x1$PC[,i.eofs]) )
      #print(dim( diag(x1$W[i.eofs]) )); print(dim( Re(t(x.m$vectors)) ))

      a.m <- t( t(x1$EOF[i.eofs,]) %*% diag(x1$W[i.eofs]) %*% Re(t(x.m$vectors)) )
      b.m <- t( t(x2$EOF[i.eofs,]) %*% diag(x2$W[i.eofs]) %*% Re(t(y.m$vectors)) )
      u.k <- x1$PC[,i.eofs] %*% Re(x.m$vectors)
      v.k <- x2$PC[,i.eofs] %*% Re(x.m$vectors)
      R <- sqrt(Re(x.m$values))
   }

    if (test) {
      # X1 = a.t%*% t(u.t)
      # X2 = b.t%*% t(v.t)
      # print(dim(a.m)); print(dim(u.k))
      X1.cca <- u.k %*% a.m
      X2.cca <- v.k %*% b.m
      dim(X1.cca) <- c(length(x1$tim),length(x1$lat),length(x1$lon))
      dim(X2.cca) <- c(length(x2$tim),length(x2$lat),length(x2$lon))
      iy <- floor(length(x1$lat)/2); ix <- floor(length(x1$lon)/2)
      jy <- floor(length(x2$lat)/2); jx <- floor(length(x2$lon)/2)
      y.1 <- EOF2field(x1,anomalies=TRUE)$dat[,iy,ix]
      y.2 <- EOF2field(x2,anomalies=TRUE)$dat[,jy,jx]
      print(summary(y.1))

      par(mfcol=c(2,1))
      plot(y.1,main="Test: CCA reconstruction",type="l",lwd=3,col="grey60")
      lines(X1.cca[,iy,ix],col="red",lty=2)

      plot(y.2,main="Test: CCA reconstruction",type="l",lwd=3,col="grey60")
      lines(X2.cca[,jy,jx],col="blue",lty=2)
      newFig()
     }
    dim(a.m) <- d.1; dim(b.m) <- d.2
 } #endif (eof)

  sub <- paste("r=",round(R[1],2),sub)
  cca <- list(a.m = a.m, b.m =b.m, u.k= u.k, v.k = v.k, r=R,
              x1=x1,x2=x2, main=main, sub=sub, i1=i1, i2=i2)
  if (test) {
    cca$C.11 <- C.11 
    cca$C.22 <- C.22 
    cca$C.12 <- C.12
    if (SVD) {cca$l.k <- l.k; cca$r.k <- r.k; cca$M <- M} 
    if (!SVD) {cca$x.m <- x.m; cca$y.m <- y.m; cca$M.1 <- M.1; cca$M.2 <- M.2} 
    cca$test.diag <- X1.cca
  }
  class(cca) <- c("CCA", class(x1)[2])

  if (round(R[1],2) != round(cor(u.k[i1,1],v.k[i2,1]),2)) {
    print("WARNING: The correlations are not internally consistent!")
    print(paste("CCA: leading canonical correlation=", round(R[1],2),
                " actual correlation=",round(cor(u.k[i1,1],v.k[i2,1]),2)))
  }

  rm(a.m, b.m, main, R, u.k, v.k, x1, x2, sub)  
  if (plot) plotCCA(cca)
  invisible(cca)
}

plotCCA <- function(cca) {
    attach(cca)
    #print(dim(t(a.m[1,,]))); print(c(length(x1$lon),length(x1$lat)))
    image(x1$lon,x1$lat,t(a.m[1,,]),col = cm.colors(21),
        main=main, sub=sub,xlim=range(c(x1$lon,x2$lon)),
        ylim=range(c(x1$lat,x2$lat)))    
    addland()
    contour(x1$lon,x1$lat,t(a.m[1,,]),lwd=2,col="darkblue",add=TRUE)
    contour(x2$lon,x2$lat,t(b.m[1,,]),lwd=1,add=TRUE,col="darkred")
    legend(min(c(x1$lon,x2$lon)),max(c(x1$lat,x2$lat)),c(x1$v.name,x2$v.name),
           col=c("darkblue","darkred"),
           lwd=c(1,2),bg="grey95")

    newFig()
    plot(x1$yy+(x1$mm-1)/12+x1$dd/365.25,u.k[,1],type="l",col="blue",
         main=paste("CCA: correlation=",round(cor(u.k[i1,1],v.k[i2,1]),2)),
         xlab="Time",ylab="Canonical variates",sub=sub)
    lines(x2$yy+(x2$mm-1)/12+x2$dd/365.25,v.k[,1],type="l",col="red",lty=2)
    legend(min(x1$yy),max(u.k[,1]),c(x1$v.name,x2$v.name),col=c("blue","red"),
           lwd=1,lty=c(1,2),bg="grey95")

    newFig()
    plot(r,pch=20,cex=1.7,main="Canonical correlations",sub=sub,ylim=c(0,1.2))
    grid()
    detach(cca)
}


MVR <- function(x,y,plot=TRUE,main="Multivariate regression",sub="",test=FALSE,i.eofs=1:8,LINPACK=TRUE,SVD=TRUE) {
# y = x %psi + %xi

  ix <- is.element(x$yy*10000 + x$mm*100 + x$dd, y$yy*10000 + y$mm*100 + y$dd)
  iy <- is.element(y$yy*10000 + y$mm*100 + y$dd, x$yy*10000 + x$mm*100 + x$dd)

  if ( (class(x)[1]=="field") & (class(y)[1]=="field") ) {
    X <- x$dat[ix,,]; d.x <- dim(X); dim(X) <- c(d.x[1],d.x[2]*d.x[3])
    Y <- y$dat[iy,,]; d.y <- dim(Y); dim(Y) <- c(d.y[1],d.y[2]*d.y[3])
        
    if (SVD) {
      if (LINPACK) UWV <-svd(X) else 
                   UWV <-La.svd(X)
      V <- UWV$v; D <- UWV$d
    #print(c(dim(t(V)),NA,dim(diag(D^2)),NA,dim(V),NA,dim(t(X)),NA,dim(Y)))
      psi <- chol2inv(t(V) %*% diag(D^2) %*% V) %*% t(X) %*% Y
    } else {
      psi <- solve(t(X) %*% X) %*% t(X) %*% Y # close to singular
    }
    Yhat <- X %*% psi
 #   Y <- t(Y); Yhat <- t(Yhat); dim(Y) <- ; dim(Yhat) <- c(126,11,31)
    mvr  <- y
    dim(Yhat) <- c(d.y[1],d.y[2],d.y[3])
    mvr$dat <- Yhat
    mvr$psi <- psi
  } else if ( (class(x)[1]=="eof") & (class(y)[1]=="eof") ) {
    d.x <- dim(x$EOF); if (length(d.x)>2) dim(x$EOF) <- c(d.x[1],d.x[2]*d.x[3])
    X <- x$PC[ix,]; Y <- y$PC[iy,]
    psi <- solve(t(X) %*% X) %*% t(X) %*% Y
    Yhat <- X %*% psi
    psi.map <- t(psi %*% diag(x$W) %*% x$EOF)
    dim(psi.map) <- c(x$size[2],x$size[3],d.x[1])

    map.x <- t(diag(x$W) %*% x$EOF)
    dim(map.x) <- c(x$size[2],x$size[3],d.x[1])
    mvr  <- y
    mvr$PC <- Yhat
    mvr$psi <- psi
    Psi <- list(lon=x$lon,lat=x$lat,map=t(psi.map[,,1]),v.name=paste('projected',y$v.name),
                tim=NULL,date=NULL,attributes=y$attributes)
    class(Psi) <- "map"
    attr(Psi,"descr") <- "MVR regression weights"
    mvr$psi <- Psi
    if (test) {
      iy <- floor(length(y$lat)/2); ix <- floor(length(y$lon)/2)
      y.1 <- EOF2field(x,anomalies=TRUE)$dat[,iy,ix]
      print(summary(c(Yhat)))
      print(class(mvr))      
      if ( (class(mvr)[1]=="eof") ) {
         print("EOF2field....")
         y.1 <- EOF2field(y,anomalies=TRUE)$dat[,iy,ix]
         y.2 <- EOF2field(mvr,anomalies=TRUE)$dat[,iy,ix]
      }
      plot(y.1,main="Test: MVR reconstruction",type="l",lwd=3,col="grey60")
      lines(y.2,col="red",lty=2)
      newFig()
      print(c(dim(Psi$map),NA,length(Psi$lon),length(Psi$lat)))
      map(Psi)

    mar.orig <- (par.orig <- par(c("mar","las","mfrow")))$mar
    on.exit(par(par.orig))

    w <- (3 + mar.orig[2]) * par('csi') * 2.54
    layout(matrix(c(2, 1), nc=2), widths=c(1, lcm(w)))
    
    par(las = 1)
    mar <- mar.orig
    mar[4] <- 1
    par(mar=mar)

    contour(x$lon,x$lat,t(map.x[,,1]),col="grey60",lwd=2,lty=2,add=TRUE)
    }      
  }
  class(mvr) <- c(class(y),"MVR")
  invisible(mvr)
}

POP <- function(x,plot=TRUE,main="POP analysis",sub="",
                test=FALSE,i.eofs=1:8,LINPACK=TRUE,mode=1) {
# After von Storch & Zwiers (1999), Statistical Analysis in Climate Research, p. 338.
  n <- length(x$tim)
   if ( (class(x)[1]=="field") ) {
    dims <- dim(X$dat)
    dim(X$dat) <- c(dims[1],dims[2]*dims[3])
    X1 <- x$dat[2:n,]; X2 <- x$dat[1:(n-1),]
    C12 <- cov(X1,X2); C11<- cov(X1,X1)
    A <- C12 %*% solve(C11)
    e <- eigen(A)
    maps <- abs(e$vectors)   
  } else if ( (class(x)[1]=="eof") ) {
    X1 <- x$PC[2:n,i.eofs]; X2 <- x$PC[1:(n-1),i.eofs]
    C12 <- cov(X1,X2); C11<- cov(X1,X1)
#    A <- C12 %*% solve(C11)
    A <- C12 / diag(C11)
    e <- eigen(A)
    # TEST e$vectors <- diag(rep(1,length(i.eofs))); print(e$vectors)
    maps <- t(abs(e$vectors) %*% x$EOF[i.eofs,])
    pop.Im <-  t(Im(e$vectors) %*% x$EOF[i.eofs,])
    pop.Re <-  t(Re(e$vectors) %*% x$EOF[i.eofs,])
    #print(dim(maps))
    dim(maps) <- c(length(x$lat),length(x$lon),length(e$values))   
    dim(pop.Im) <- c(length(x$lat),length(x$lon),length(e$values))   
    dim(pop.Re) <- c(length(x$lat),length(x$lon),length(e$values))   
 }
  pop <- e
  pop$maps <- maps
  pop$pop.Im <- pop.Im
  pop$pop.Re <- pop.Re
  pop$lon <- x$lon; pop$lat <- x$lat
  class(pop) <- c("POP",class(x))
  if (plot) plotPOP(pop)
  invisible(pop)
}

plotPOP <- function(pop,mode=1,main="POP analysis",sub="") {
     if ( (class(pop)[1]!="POP") ) stop("Need a 'POP' object")

     newFig()
     plot(pop$values,pch=20,col="grey50",main=main,sub=sub)
     points(pop$values)
     grid()

     if (sub=="") sub <- paste("Real=",round(Re(pop$values[mode]),2),
                               "Imaginary=",round(Im(pop$values[mode]),2))
     main <- paste(main," Mode=",mode)
     newFig()
     image(pop$lon, pop$lat,t(pop$maps[,,mode]),main=main,sub=sub)
     grid()
     addland()
     contour(pop$lon, pop$lat,t(pop$pop.Re[,,mode]),lwd=2,col="grey40",add=TRUE)
     contour(pop$lon, pop$lat,t(pop$pop.Im[,,mode]),lty=2,add=TRUE)
}


SSA <- function(x,m,plot=TRUE,main="SSA analysis",sub="",
                LINPACK=TRUE,param="t2m",anom=TRUE) {
# After von Storch & Zwiers (1999), Statistical Analysis in Climate Research, p. 312.  

  if (class(x)[1] != "station") stop('SSA: need a station object')
  x.mean <- 0

  if (class(x)[2] == "monthly.station.record") param <- "val"

  if (anom) x <- anomaly.station(x) else {
    expr <- paste("x$",param,sep="")
    x.mean <- mean(eval(parse(text=expr)),na.rm=TRUE)
    #print(x.mean)
  }

  if (class(x)[2] == "monthly.station.record") {
    nt <- length(x$yy)*12 
    x$val <- c(t(x$val))
  } else if (class(x)[2] == "daily.station.record") {
    nt <- length(x$yy)
  }

  Nm <- nt - m + 1
  X <- matrix(rep(NA,Nm*m),Nm,m)
  #print(dim(X))
  for (i in 1:m) {
    expr <- paste("x$",param,"[",i,":(nt-",m-i,")] - x.mean",sep="")
    #print(expr); #print(length(eval(parse(text=expr))))
    X[,i] <- eval(parse(text=expr)) 
  }

  #image(1:dim(X)[1],1:dim(X)[2],X)
   
  if (LINPACK) ssa <- svd(X) else
               ssa <- La.svd(X)

  if (sub=="") sub <- paste("Window width=",m)

  ssa$m<- m; ssa$Nm <- Nm; ssa$nt <- nt
  ssa$anom <- anom
  ssa$param <- x$param
  ssa$station <- x
  class(ssa) <- c("SSA",class(x))
  if (plot) plotSSA(ssa)
  invisible(ssa)
}

lagStation <- function(obs,lag=1) {
  dims <- dim(t(obs$val))
  ts <- c(t(obs$val)); n <- length(ts)
  if (lag<0) for (i in 1:abs(lag)) ts <- c(ts[-1],NA) else
             for (i in 1:abs(lag)) ts <- c(NA,ts[-n])
  dim(ts) <- dims
  obs$val <- t(ts)
  class(obs) <- c(class(obs),paste("lagged",lag))
  invisible(obs)
}

plotSSA <- function(ssa,main="SSA analysis",sub="")  {
    if ( (class(ssa)[1]!="SSA") ) stop("Need an 'SSA' object")
    nt <- ssa$nt
    newFig()
    plot(ssa$d,main=main,sub=sub,ylab="Singular value",pch=20,col="grey50")
    points(ssa$d)
    grid()

    newFig()
    par(mfcol=c(3,1))
    plot(ssa$v[,1],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA vector: mode 1",lwd=3,col="grey70")
    grid()
    plot(ssa$v[,2],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA vector: mode 1",lwd=3,col="grey70")
    grid()
    plot(ssa$v[,3],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA vector: mode 1",lwd=3,col="grey70")
    grid()


    newFig()
    par(mfcol=c(3,1))
    if (class(ssa)[3] == "monthly.station.record") {
      yy <- sort(rep(ssa$x$yy,12)); yy <- yy[1:ssa$Nm]
      mm <- rep(1:12,nt); mm <- mm[1:ssa$Nm]
      #print(dim(ssa$v)); print(dim(ssa$u)); print(length(yy)); print(length(mm));
      plot(yy + (mm-0.5)/12, ssa$u[,1],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
      grid()
      plot(yy + (mm-0.5)/12, ssa$u[,2],type="l",main=main,syb=sub,
           xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
      grid()
      plot(yy + (mm-0.5)/12, ssa$u[,3],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
      grid()
    } else if (class(ssa)[3] == "daily.station.record") {
      plot(ssa$x$yy[1:ssa$Nm] + ssa$x$mm[1:ssa$Nm]/12 + ssa$x$dd[1:ssa$Nm]/365, ssa$u[,1],
           type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
      grid()
      plot(ssa$x$yy[1:ssa$Nm] + ssa$x$mm[1:ssa$Nm]/12 + ssa$x$dd[1:ssa$Nm]/365, ssa$u[,2],
           type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
      grid()
      plot(ssa$x$yy[1:ssa$Nm] + ssa$x$mm[1:ssa$Nm]/12 + ssa$x$dd[1:ssa$Nm]/365, ssa$u[,3],
           type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
      grid()
    }
  }

testCCA <- function(method="CCA",reconstr=FALSE,mode=1,test=TRUE,LINPACK=TRUE,SVD=TRUE,n.pc=4,synthetic=TRUE) {
  print("version 0.1:")
  data(eof.slp,envir=environment())
  print(dim(eof.slp$EOF)); print(dim(eof.slp$PC)); print(length(eof.slp$W))
  eof.slp$EOF[!is.finite(eof.slp$EOF)] <- 0
  print(summary(c(eof.slp$EOF)))
  print(summary(c(eof.slp$PC)))
  print(summary(c(eof.slp$W)))
  eof.slp$clim <- eof.slp$EOF[1,]*0
  dim(eof.slp$clim) <- c(73,144)
  if (synthetic) {
    nt <- 200
    eof.slp$tim <- 1:nt; eof.slp$yy <- 2000 + floor((1:nt)/360)
    eof.slp$mm <- mod(floor((1:nt)/30),12)+1; eof.slp$dd <- mod(1:nt,30)+1 
    eof.slp$size[1,1] <- nt
    #print(rbind(eof.slp$tim,eof.slp$yy,eof.slp$mm,eof.slp$dd))
    eof.slp$PC <- matrix(rep(0,n.pc*nt),nt,n.pc)
    eof.slp$id.t <- rep("test",nt)
  } else nt <- length(eof.slp$tim)
  eof1 <- eof.slp; eof2 <- eof.slp; rm(eof.slp)
  eof1$PC <- eof1$PC[,1:n.pc]; eof1$EOF <- eof1$EOF[1:n.pc,]; eof1$W <- eof1$W[1:n.pc]
  eof2$PC <- eof2$PC[,1:n.pc]; eof2$EOF <- eof2$EOF[1:n.pc,]; eof2$W <- eof2$W[1:n.pc]
  eof1$dW <- eof1$dW[1:n.pc]; eof2$dW <- eof2$dW[1:n.pc]
  eof1$var.eof <- eof1$var.eof[1:n.pc]; eof2$var.eof <- eof2$var.eof[1:n.pc]

  print(dim(eof1$EOF)); print(dim(eof1$PC)); print(length(eof1$W))
  y.1 <- EOF2field(eof1,anomalies=TRUE)$dat[,1,1]
  print(summary(y.1))
  
  if (synthetic) {
    modes <- 1:n.pc; modes <- modes[-mode]
    print("signal: sin")
    eof1$PC[,mode] <- 50*sin(seq(-12*pi,12*pi,length=nt))
    eof2$PC[,mode] <- 50*sin(seq(-12*pi,12*pi,length=nt))
    print("noise: rnorm")
    for (i in modes) {
      eof1$PC[,i] <- rnorm(nt)
      eof2$PC[,i] <- rnorm(nt)
    }
  }
  if (reconstr) {
    print("Reconstruct fields...")
    x1 <- EOF2field(eof1)
    x2 <- EOF2field(eof2)
    print(class(x1))
    print("Run test...")
    print(paste(method,"(x1,x2,test=",test,
                 ",LINPACK=",LINPACK,",SVD=",SVD,")",sep=""))
    cca.test <- eval(parse(text=paste(method,"(x1,x2,test=",test,
                 ",LINPACK=",LINPACK,",SVD=",SVD,")",sep="")))
  } else {
    print(paste(method,"(eof1,eof2,test=",test,
                 ",LINPACK=",LINPACK,",SVD=",SVD,")",sep=""))
    cca.test <- eval(parse(text=paste(method,"(eof1,eof2,test=",test,
                 ",LINPACK=",LINPACK,",SVD=",SVD,")",sep="")))
  }
  invisible(cca.test)
}

Psi <- function(cca) {
  G <- cca$a.m; d1 <- dim(G); dim(G) <- c(d1[1],d1[2]*d1[3])
  H <- cca$b.m; d2 <- dim(H); dim(H) <- c(d2[1],d2[2]*d2[3])
  M <- diag(cca$r)
  V <- cca$u.k
  U <- cca$v.k
  G <- t(G); H <- t(H)
 
  # print(dim(G)); print(dim(M)); print(dim(H))
  if (class(cca)[1] =="CCA") Psi <- G %*% M %*% solve(t(H) %*% H) %*% t(H)
  #if (class(cca)[1] =="SVD") Psi <- G %*% M %*% solve(Cxx) %*% t(H)
  class(Psi) <- paste(class(cca)[1],"model",sep=".")
  attr(Psi,"dims") <- d1
  attr(Psi,"lon") <- cca$x1$lon
  attr(Psi,"lat") <- cca$x1$lat
  Psi
}
  
predictCCA <- function(Psi,X) {
  if ( (class(X)[1]!="eof") & (class(X)[1]!="field")) stop('Need a field or EOF object!')
  type <- class(X)
  if (type[1]=="eof") field <- EOF2field(X)
  X <- field$dat
  d <- dim(X); dim(X) <- c(d[1],d[2]*d[3])
  X <- t(X)
  #print(dim(Psi)); print(dim(X)); print(d)
  Y.hat <-  Psi %*% X
  field$dat <- t(Y.hat)
  #print(dim(field$dat))
  d1 <- attr(Psi,"dims")
  dim(field$dat) <- c(d[1],d1[2],d1[3])
  field$lon <- attr(Psi,"lon"); nx <- length(field$lon)
  field$lat <- attr(Psi,"lat"); ny <- length(field$lat)
  field$id.x <- rep("CCA",nx*ny)
  field$id.lon <- rep("CCA",nx)
  field$id.lat <- rep("CCA",ny)
  field$id.t <- rep("CCA",d[1])
  #print("HERE")
  if (type[1]=="eof") result <- EOF(field) else result <- field
  result
}


check.repeat <- function(x) {

  if (max(as.numeric(table(x$yy)))==1) return(x)
  
  yy <- seq(min(x$yy),max(x$yy),by=1); ny <- length(yy)
  val <- matrix(rep(NA,ny*12),ny,12)
  for (it in 1:ny) {
    ii <- (1:ny)[is.element(yy[it],x$yy)]
    #print(yy[it]); print(ii); print(ii[1]); print(dim(val)); print(dim(x$val))
    val[it,] <- x$val[ii[1],]
  }
  x$val <- val
  x$yy <- yy
  invisible(x)
}

stations2field <- function(data.set=c("narp"),ele=101,obj.type="monthly.field.object",
                           plot=TRUE,silent=FALSE,intrp.method="interpp",
                           interpolation.option="simple") {
  first.obs <- TRUE
  if (!silent) print("stations2field")

  if (is.character(data.set)) {
    for (i.data in 1:length(data.set)) {
      method <- paste("get",data.set[i.data],sep="")
      if (!silent) print(paste(method,"()",sep=""))
      stations<- eval(parse(text=paste(method,"()",sep="")))
      if (is.list(stations)) {
        if(!is.null(stations$number)) stations <-  as.numeric(stations$number) else 
        if(!is.null(stations$stnr)) stations <-  as.numeric(stations$stnr) else
        if(!is.null(stations$station)) stations <-  as.numeric(stations$station) else
        if(!is.null(stations$name)) stations <-  stations$name else
        if(!is.null(stations$location)) stations <-  stations$location
      }
      if (!silent) print(stations)
      ns <- length(stations)
      is <- 1
      while (is <= ns) {
        if (!silent) print(stations[is])
        if (is.character(stations[is])) obs <-
               eval(parse(text=paste(method,"('",stations[is],"',ele=",ele,")",sep=""))) else
        if (is.numeric(stations[is])) obs <-
               eval(parse(text=paste(method,"(",stations[is],",ele=",ele,")",sep="")))
        if (length(obs$yy) > 0)  {
          obs <- check.repeat(obs)
          if (first.obs) {
            yy.int <- range(obs$yy)
            yy <- sort(rep(yy.int[1]:yy.int[2],12))
            mm <- rep(1:12,length(yy.int[1]:yy.int[2]))
            nt <- length(yy)
            Dat <- matrix(rep(NA,ns*nt),ns,nt)
            Lon <- rep(NA,ns); Lat <- Lon
            first.obs <- FALSE
          }
          yy.obs <- sort(rep(obs$yy,12))
          mm.obs <- rep(1:12,length(obs$yy))
          y <- c(t(obs$val))
          i1 <- is.element(yy*100+mm,yy.obs*100+mm.obs)
          i2 <- is.element(yy.obs*100+mm.obs,yy*100+mm)
          Dat[is,i1] <- y[i2]
          Lon[is] <- obs$lon
          Lat[is] <- obs$lat
         } else {
          ns <- ns-1
       }
       is <- is+1
      }
    }
  } else {
    print("Presumes that 'data.set' is a list of station objects")
    stations<- names(data.set)
    ns <- length(stations)
    for (is in 1:ns) {
        if (!silent) print(stations[is])
        obs <- eval(parse(text=paste("data.set$",stations[is],sep="")))
        obs <- check.repeat(obs)
        if (first.obs) {
          yy.int <- range(obs$yy)
          yy <- sort(rep(yy.int[1]:yy.int[2],12))
          mm <- rep(1:12,length(yy.int[1]:yy.int[2]))
          nt <- length(yy)
          Dat <- matrix(rep(NA,ns*nt),ns,nt)
          Lon <- rep(NA,ns); Lat <- Lon
          first.obs <- FALSE
        }
        yy.obs <- sort(rep(obs$yy,12))
        mm.obs <- rep(1:12,length(obs$yy))
        y <- c(t(obs$val))
        i1 <- is.element(yy*100+mm,yy.obs*100+mm.obs)
        i2 <- is.element(yy.obs*100+mm.obs,yy*100+mm)
        Dat[is,i1] <- y[i2]
        Lon[is] <- obs$lon
        Lat[is] <- obs$lat
    }
  }

  i.val <- is.finite(Lon)
  Dat <- Dat[i.val,]
  Lon <- Lon[i.val]
  Lat <- Lat[i.val]
  i.val <- is.finite(colMeans(Dat))
  Dat <- Dat[,i.val]
  yy <- yy[i.val]; nt <- length(yy)
  mm <- mm[i.val]
  x.centre <- mean(Lon)
  y.centre <- mean(Lat)
  xy <- COn0E65N(Lon, Lat,lon.0=x.centre,lat.0=y.centre)

  lat <- seq(min(Lat),max(Lat),by=0.25); ny <-length(lat)
  lon <- seq(min(Lon),max(Lon),by=0.50); nx <-length(lon)
  lon.xy <- rep(lon,length(lat))
  lat.xy <- sort(rep(lat,length(lon)))
  XY <- COn0E65N(lon.xy, lat.xy,lon.0=x.centre,lat.0=y.centre)
  
  dat <- matrix(rep(NA,nt*nx*ny),nt,ny*nx); dim(dat) <- c(nt,ny,nx)
  if (!silent) print(c(dim(dat),NA,length(XY$x),NA,length(lon),length(lat),NA,
                         nx,ny,nx*ny,length(lon.xy)))

  interpolation.approach <- switch(interpolation.option,
     "simple"=paste("interp(Lon[ok],Lat[ok],Dat[ok,i],lon,lat,duplicate='mean')",sep=""),
     "distance"=paste(intrp.method,"(xy$x[ok],xy$y[ok],Dat[ok,i],XY$x,XY$y,duplicate='mean')",sep=""),
     "test"=paste(intrp.method,"(Lon[ok],Lat[ok],Dat[ok,i],lon.xy,lat.xy,duplicate='mean')",sep=""))

  for (i in 1:nt) {
    ok <- is.finite(Dat[,i])
    map <- eval(parse(text= interpolation.approach))
    if (is.list(map)) map <- map$z
    dim(map) <- c(nx,ny)
    if(plot) {
      plot(range(lon),range(lat),type="n",main=paste(i,nt),
           sub=intrp.method)
      addland()
      image(lon,lat,map,add=TRUE,col = cm.colors(21))
      contour(lon,lat,map,add=TRUE)
      points(Lon,Lat)
      text(Lon,Lat,round(Dat[,i],1))
    }
    dat[i,,] <- t(map)
  }

  v.nam <- switch(as.character(ele),"101"="T2m","601"="Precip")
  data.sets <- ""
  for (i.data in 1:length(data.set)) data.sets <- paste(data.sets,data.set[i],sep="+")
  id.x <- rep(v.nam,nx*ny); dim(id.x) <- c(ny,nx)
  id.t <- rep(v.nam,nt)
  tim <- 1:nt
  attr(tim,"unit") <- "month"
  attr(tim,"time_origin") <- paste(15,mm[1],yy[1],sep="-")
  dat.att <- list(time.unit="month", time.origin=paste(15,mm[1],yy[1],sep="-"),
                  unit=obs$unit,
                  long.name=paste("Field object created from",data.set,obs$obs.name),
                  filename=NULL,scale.factor=1,add.offset=0,miss=NA,
                  daysayear=365.25)
  field  <- list(dat=dat,lon=lon,lat=lat,tim=tim,lev=NULL,
                 v.name=v.nam,id.x=id.x,id.t=id.t,
                 yy=yy,mm=mm,dd=rep(15,length(yy)),n.fld=1,
                 id.lon=rep(v.nam,nx),id.lat=rep(v.nam,ny),
                 attributes=dat.att,
                 filename=NULL,Lon.src=Lon,Lat.src=Lat)
  class(field) <- c("field",obj.type)
  invisible(field)
}
