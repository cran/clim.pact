# Computes approximate cartesian coordinates in km
# centered on 0E 65N.
# R.E. Benestad, DNMI, 04.01.2001
#
COn0E65N <- function(lon, lat, lat.0=65,lon.0=0) {
  a<-6357 # km
  if (length(lon) != length(lat)) {
    stop("lon and lat must have same length")
  }
  COn0E65N<-list(y=a * pi*(lat-lat.0)/180,
                 x=a*pi*(lon-lon.0)*cos(pi*lat/180)/180)
  COn0E65N
  }
