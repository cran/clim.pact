load("koebenhavn.t2m.Rdata")
koebenhavn.t2m$location <- "Copenhagen"
class(koebenhavn.t2m) <- c("station","monthly.station.record")
dat <- koebenhavn.t2m
