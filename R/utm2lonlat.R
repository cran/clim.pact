# R.E. Benestad, met.no, Oslo, Norway 18.08.2003
# rasmus.benestad@met.no
#------------------------------------------------------------------------

utm2lonlat <- function(easting,northing,zone,band=NULL,a=6.378e06) {

  print("S I M P L E and A P P R O X I M A T E UTM-to-Lon/Lat conversion")
  print("Assuming a sphere")
  
  zones <- seq(-177,177,by=6)
  bands <- seq(-80,80,by=8)
  bands[21] <- 84
  bandnms <- c("C","D","E","F","G","H","J","K","L","M","N","P",
               "Q","R","S","T","U","V","W","X")
  cent.mer <- min(zones[zone >= zones])
  theta <- 180*(easting-500000)/(a*pi) + cent.mer
  phi <- 180*northing/(a*pi)
  if (!is.null(band)) {
    ii <- is.element(bandnms,uppercase(band))
    if (bands[ii] < 0) {
      phi <- 180*(northing-10000000)/(a*pi)
      phi <- -phi
    }
  }
  lonlat=list(lon=theta,lat=phi,easting=easting,
              northing=northing,zone=zone,band=band)
  class(lonlat) <- "LonLat.coordinate"
  invisible(lonlat)
}
