load("helsinki.t2m.Rdata")
helsinki.t2m$location <- "Helsinki"
class(helsinki.t2m) <- c("station","monthly.station.record")
dat <- helsinki.t2m
