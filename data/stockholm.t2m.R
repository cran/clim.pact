load("stockholm.t2m.Rdata")
stockholm.t2m$location <- "Stockholm"
class(stockholm.t2m) <- c("station","monthly.station.record")
dat <- stockholm.t2m
