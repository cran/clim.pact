distAB <- function(lon,lat,lons,lats,a=6.378e06) {
  theta <- pi*lon/180
  phi <- pi*lat/180
  dist <- rep(NA,length(lons))
  r1 <- c(cos(phi)*cos(theta),
          sin(phi),
          cos(phi)*sin(theta))
  for (i in 1:length(lons)) {
    theta <- pi*lons[i]/180
    phi <- pi*lats[i]/180
    r2 <- c(cos(phi)*cos(theta),
            sin(phi),
            cos(phi)*sin(theta))
#    angle <- acos( sum(r1*r2)/(sqrt(sum(r1*r1)) * sqrt(sum(r2*r2))) )
    angle <- acos( sum(r1*r2) )
    dist[i] <- a* angle
  }
  dist
}
                   
