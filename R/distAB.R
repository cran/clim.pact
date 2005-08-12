distAB <- function(lon,lat,lons,lats,a=6.378e06) {
  lons <- lons[is.finite(lons)]
  lats <- lats[is.finite(lats)]
  theta <- pi*lon/180
  phi <- pi*lat/180
  dist <- rep(NA,length(lons))
  r1 <- c(cos(phi)*cos(theta),
          sin(phi),
          cos(phi)*sin(theta))
  dim(r1) <- c(3,length(lon))
  theta <- pi*lons/180
  phi <- pi*lats/180

  r2 <- cbind(cos(phi)*cos(theta),
              sin(phi),
              cos(phi)*sin(theta))
#  angle <- acos( sum(r1*r2) )
  angle <- acos( r2 %*% r1 )
#  if (sum(!is.finite(angle))>0) {
#    ibad <- !is.finite(angle)
#    print("MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM: distAB")
#    print("Detected a not-a-finite number...")
#    print(r2[ibad,])
#    print(r1)
#    print("MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM: distAB")
#  }
  dist <- a* angle
  dist
}
                   
