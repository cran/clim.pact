# This routine computes the month, day, and year, given a Julian day.
# The algorithm is taken from Press et al. (1989), "Numerical Recipes 
# in Pascal", Cambridge, p. 13.
#
# This function removes the dependency to outdated packages 'chron' and
# 'date'.
#
# R.E. Benestad, met.no, Oslo, Norway 04.09.2003
# rasmus.benestad@met.no
#------------------------------------------------------------------------

caldat <- function(julian) {

  igreg=2299161
  julian <- trunc(julian)
  jalpha <- julian*0; ja <- julian*0
  im <-  (julian >= igreg)
  if (sum(im)>0) {
    jalpha[im] <- trunc(((julian-1867216) - 0.25)/36524.25) # Cross-over to Gregorian Calendar
                                                        # produces this correction.
    ja[im] <- julian+1+jalpha-trunc(0.25*jalpha)
  }
  im <- (julian < igreg)
  if (sum(im)>0) ja[im] <- julian[im]
  jb <- ja + 1524
  jc <- trunc(6680 + ((jb-2439870)-122.1)/365.25)
  jd <- 365*jc+trunc(0.25*jc)
  je <- trunc((jb-jd)/30.6001)
  id <- jb-jd-trunc(30.6001*je)
  mm <- je-1
  im  <- (mm > 12)
  if (sum(im)>0) mm[im] <- mm[im] - 12
  iyyy <- jc-4715
  im <-  (mm > 2)
  if (sum(im)>0) iyyy[im] <- iyyy[im] - 1
  im <-  (iyyy <= 0)
  if (sum(im)>0) iyyy <- iyyy - 1
  
  caldat <- list(month=mm,day=id,year=iyyy)
  invisible(caldat)
}
