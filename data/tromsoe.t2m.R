load("tromsoe.t2m.Rdata")
tromsoe.t2m$location <- "Tromsoe"
class(tromsoe.t2m) <- c("station","monthly.station.record")
dat <- tromsoe.t2m
