# R.E. Benestad, met.no, Oslo, Norway 18.08.2003
# rasmus.benestad@met.no
#------------------------------------------------------------------------

lonlat2utm <- function(lon,lat,a=6.378e06) {

  print("S I M P L E and A P P R O X I M A T E Lon/Lat-to-UTM conversion")
  print("Assuming a sphere")
  
  zones <- seq(-177,177,by=6)
  zonesnms <- 1:length(zones)
  bands <- seq(-80,80,by=8)
  bands[21] <- 84
  bandsnms <- c("C","D","E","F","G","H","J","K","L","M","N","P",
                "Q","R","S","T","U","V","W","X")
  zone <- zonessnms[lon >= zones][1]
  band <- bandsnms[lat >= bands][1]
  easting  <-  pi*(zones[lon >= zones][1] - lon)/180*a + 500000
  if (lat > 0) northing <- pi*lat/180*a else
               northing <- -pi*lat/180*a + 10000000
  utc <- list(zone=zone,band=band,northing=northing,easting=easting,lat=lat,lon=lon)
  class(utc) <- "UTC.coordinate"
  invisible(utc)
}
