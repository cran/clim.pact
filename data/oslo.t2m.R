load("oslo.t2m.Rdata")
oslo.t2m$location <- "Oslo"
class(oslo.t2m) <- c("station","monthly.station.record")

