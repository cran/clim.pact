load("bergen.t2m.Rdata")
bergen.t2m$location <- "Bergen"
class(bergen.t2m) <- c("station","monthly.station.record")
dat <- bergen.t2m
