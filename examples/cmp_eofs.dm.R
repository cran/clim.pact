rm(list=ls())
gc()
library(clim.pact)

reduce.res <- function(x,grid.skip=NULL,max.grids=20) {
  nx <- length(x$lon); ny <- length(x$lat)
  if (!is.null(grid.skip)) {
    x.grid.sel <- seq(1,nx,by=grid.skip)
    y.grid.sel <- seq(1,ny,by=grid.skip)
  } else if (!is.null(max.grids)) {
    x.grid.sel <- round(seq(1,nx,length=max.grids))
    y.grid.sel <- round(seq(1,ny,length=max.grids))
  }
  x$dat <- x$dat[,y.grid.sel,x.grid.sel]
  x$lon <- x$lon[x.grid.sel]
  x$lat <- x$lat[y.grid.sel]
  x$id.lon <- x$id.lon[x.grid.sel]
  x$id.lat <- x$id.lat[y.grid.sel]
  x$id.x <- x$id.x[y.grid.sel,x.grid.sel]
  invisible(x)
}

# Domains:
# x.rng <- c(3,12);  y.rng <- c(57,65)   # Southern Norway focus
# x.rng <- c(16,31); y.rng <- c(64,73)   # Northern Norway focus
# x.rng <- c(-5,12); y.rng <- c(55,65)   # Western Norway focus
# x.rng <- c(4,25);  y.rng <- c(58,70)   # The whole country focus
# x.rng <- c(12,26); y.rng <- c(70,80)   # Svalbard focus

domains <- c("D10","D6","D7","D8","D1","D2","D3,","D4","D5","D9")
#domains <- c("D10","D3")
print(domains)

for (domain in domains) {
  x.domain <- switch(domain,
                "D1"=list(x.rng=c(5,12),y.rng=c(57,66)),
                "D2"=list(x.rng=c(20,31),y.rng=c(66,72)),
                "D3"=list(x.rng=c(3,9),y.rng=c(56,63)),
                "D4"=list(x.rng=c(12,26),y.rng=c(74,80)),
                "D5"=list(x.rng=c(7,20),y.rng=c(60,70)),
                "D6"=list(x.rng=c(4,14),y.rng=c(58,64)),
                "D7"=list(x.rng=c(25,32),y.rng=c(70,75)),
                "D8"=list(x.rng=c(7,15),y.rng=c(56,62)),
                "D9"=list(x.rng=c(4,25),y.rng=c(58,70)),
                "D10"=list(x.rng=c(7,12),y.rng=c(59,62)))
x.rng <- x.domain$x.rng
y.rng <- x.domain$y.rng
                                        # Daily values:
x.1.dm<-retrieve.nc("/data1/era15/ERA-15_t2m.nc",
                    x.rng=x.rng,y.rng=y.rng)

#EOF(x.1.dm,mon=1)
#gc()  
#EOF(x.1.dm,mon=2)
#gc()  
#EOF(x.1.dm,mon=3)
#gc()  
#EOF(x.1.dm,mon=4)
#gc(verbose=TRUE)  
#EOF(x.1.dm,mon=NULL)  
X.1.dm<-retrieve.nc("/data1/hirham/T2M_198001-199912.nc",
                    x.rng=x.rng,y.rng=y.rng)
Y.1.dm<-retrieve.nc("/data1/hirham/T2M_203001-204912.nc",
                    x.rng=x.rng,y.rng=y.rng)
Y.1.dm$yy <- Y.1.dm$yy + 50

xX.1.dm <- catFields(X.1.dm,Y.1.dm,demean=FALSE)
xX.1.dm <- catFields(x.1.dm,xX.1.dm)
xX.1.dm <- reduce.res(xX.1.dm)  
#EOF(xX.1.dm,mon=1)
#gc()  
#EOF(xX.1.dm,mon=2)
#gc()  
#EOF(xX.1.dm,mon=3)
#gc()  
#EOF(xX.1.dm,mon=4)
rm(X.1.dm,Y.1.dm,x.1.dm)
print("mixed- EOF!")
print(dim(xX.1.dm))
gc(verbose=TRUE)  
EOF(xX.1.dm,mon=NULL)

x.2.dm<-retrieve.nc("/data1/era15/ERA-15_slp.nc",
                    x.rng=x.rng,y.rng=y.rng)
#EOF(x.2.dm,mon=1)
#EOF(x.2.dm,mon=2)
#EOF(x.2.dm,mon=3)
#EOF(x.2.dm,mon=4)
#gc()
  
X.2.dm<-retrieve.nc("/data1/hirham/PSL_198001-199912.nc",
                    x.rng=x.rng,y.rng=y.rng)
Y.2.dm<-retrieve.nc("/data1/hirham/PSL_203001-204912.nc",
                    x.rng=x.rng,y.rng=y.rng)
Y.2.dm$yy <- Y.2.dm$yy + 50
xX.2.dm <- catFields(X.2.dm,Y.2.dm,demean=FALSE)
xX.2.dm <- catFields(x.2.dm,xX.2.dm)
xX.2.dm <- reduce.res(xX.2.dm)  
EOF(xX.2.dm,mon=1)
EOF(xX.2.dm,mon=2)
EOF(xX.2.dm,mon=3)
EOF(xX.2.dm,mon=4)
EOF(xX.2.dm,mon=NULL)
gc()

xX.1.dm <- reduce.res(xX.1.dm)
xX.2.dm <- reduce.res(xX.2.dm)
xX.dm <- mixFields(xX.1.dm,xX.2.dm,mon=c(12,1,2))
EOF(xX.dm,mon=1)
gc()
  
xX.dm <- mixFields(xX.1.dm,xX.2.dm,mon=c(3,4,5))
EOF(xX.dm,mon=2)
gc()
  
xX.dm <- mixFields(xX.1.dm,xX.2.dm,mon=c(6,7,8))
EOF(xX.dm,mon=3)
gc()
  
xX.dm <- mixFields(xX.1.dm,xX.2.dm,mon=c(9,10,11))
EOF(xX.dm,mon=4)
gc()
  
xX.dm <- mixFields(xX.1.dm,xX.2.dm,mon=1:12)
EOF(xX.dm,mon=NULL)
gc()
  
rm(xX.dm,x.1.dm,x.2.dm,X.1.dm,X.2.dm,Y.1.dm,Y.2.dm)  
}

q(save="no")
