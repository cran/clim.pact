# This routine computes the Julian day given a month, day, and year.
# The algorithm is taken from Press et al. (1989), "Numerical Recipes 
# in Pascal", Cambridge, p. 10.
#
# This function removes the dependency to outdated packages 'chron' and
# 'date'.
#
# R.E. Benestad, met.no, Oslo, Norway 04.09.2003
# rasmus.benestad@met.no
# Bug correction. 04.02.2005: 'jy[im] <- iyyy' -> 'jy[im] <- iyyy[im]'
# 'jm[im]' <- 'mm+1' -> 'jm[im] <- mm[im]+1', 
# 'jy[im] <- iyyy-1' -> 'jy[im] <- iyyy[im]-1'
# 'jm[im] <- mm+13' -> 'jm[im] <- mm[im]+13'
# Previous warnings: 'number of items to replace is not a multiple of replacement length'
#------------------------------------------------------------------------

julday <- function(mm,id,iyyy) {
  igreg <- 588829
  mm <- trunc(mm)
  id <- trunc(id)
  iyyy <- trunc(iyyy)
  im <-  (iyyy == 0)
  if (sum(im)>0) return("There is no year zero!")
  if ((length(mm) != length(id)) | (length(mm) != length(iyyy)) |
      (length(iyyy) != length(id))) return("The vectors must have same length!")
  im <-  (iyyy < 0)
  if (sum(im)>0) iyyy[im] <- iyyy[im]+1
  jy <- mm*0; jm <- mm*0; ja <- mm*0; 
  im <-  (mm > 2)
  if (sum(im)>0) {
    jy[im] <- iyyy[im]
    jm[im] <- mm[im]+1
  }
  im <-  (mm <= 2)
  if (sum(im)>0) {
    jy[im] <- iyyy[im]-1
    jm[im] <- mm[im]+13
  }
  jul <- trunc(365.25*jy) + trunc(30.6001*jm) + id + 1720995
  im <- (id+31*(mm+12*iyyy)>= igreg)
  if (sum(im)>0) {
    ja[im] <- trunc(0.01*jy)
    jul[im] <- jul+2-ja[im]+trunc(0.25*ja[im])
  }
  julday <- jul
  invisible(julday)
}
