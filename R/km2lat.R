# Computes approximate cartesian coordinates in km
# centered on 0E 65N.
# R.E. Benestad, DNMI, 04.01.2001
#
km2lat <- function(x, y, x.centre=0, y.centre=65) {
  a<-6357 # km
  km2lat<-180/pi * y /a + y.centre
  km2lat
  }
