rm(list=ls())

library(clim.pact)
x.1 <- retrieve.nc("/home/kareb/data/ncep/ncep_t2m.nc",
                   x.rng=c(-60,40),y.rng=c(50,75))
x.2 <- retrieve.nc("/home/kareb/data/ncep/ncep_slp.nc",
                   x.rng=c(-60,40),y.rng=c(50,75))
print(x.1$v.name)

print("Read GCM predictor data.")
X.1 <- retrieve.nc("data/mpi-gsdio_t2m.nc",
                   x.rng=c(-60,40),y.rng=c(50,75))
X.2 <- retrieve.nc("data/mpi-gsdio_slp.nc",
                   x.rng=c(-60,40),y.rng=c(50,75))
print(X.1$v.name)
print("Cat fields.")
xX.1 <- cat.fields(x.1,X.1,interval.1=c(1958,1998),interval.2=c(1958,2050))
xX.2 <- cat.fields(x.2,X.2,interval.1=c(1958,1998),interval.2=c(1958,2050))
xX <- mix.fields(xX.1,xX.2,mon=1,
                 interval=c(1900,2050))
print("EOF")
eof.c <- eof(xX.1,mon=1)
eof.mc <- eof(xX,mon=1)
eof.mc2 <- eof(xX,mon=1,lon=c(0,40),lat=c(55,69))
save(file="eof.c.Rdata",eof.c)
save(file="eof.mc.Rdata",eof.mc)
save(file="eof.mc2.Rdata",eof.mc2)





