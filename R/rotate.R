rotate <- function(lons,lats,lon.0=NULL,lat.0=NULL) {

  if (is.null(lon.0)) lon.0 <- mean(lons)
  if (is.null(lat.0)) lat.0 <- mean(lats)
  
  theta0 <- pi*lon.0/180
  phi0 <- pi*lat.0/180
  r0 <- c(cos(phi0)*cos(theta0),
          sin(phi0),
          cos(phi0)*sin(theta0))

  theta <- rep(NA,length(lons)); phi <- theta
  for (i in 1:length(lons)) {
    thetaA <- pi*lons[i]/180
    phiA <- pi*lats[i]/180
    r <- c(cos(phiA)*cos(thetaA),
           sin(phiA),
           cos(phiA)*sin(thetaA))
    thetaA <- pi*lon.0/180
    r1 <- c(cos(phiA)*cos(thetaA),
            sin(phiA),
            cos(phiA)*sin(thetaA))
    a <- r - r0
    b <- r1 - r0
    phi[i] <- acos( sum(r*r0) )
    theta[i] <- acos(sum(a*b) / (sqrt(sum(a*a)) * sqrt(sum(b*b))) )
 }
  
#  thetaA <- pi*thetaA/180
#  dtheta <- pi*dtheta/180
#  phiA <- pi*phiA/180
#  dphi <- pi*dphi/180
#    
#  phi <- rep(NA,length(thetaA))
#  theta <- rep(NA,length(thetaA))
#
#  for (i in 1:length(lons)) {
#    tA <- thetaA[i]
#    R <- rbind(
#     c(cos(dphi)*cos(dtheta),-cos(dphi)*sin(dtheta),sin(dphi)*cos(tA+dtheta)),
#     c(cos(dphi)*sin(dtheta),cos(dphi)*cos(dtheta),sin(dphi)*sin(tA+dtheta)),
#     c(-sin(dphi-tA),cos(dphi-tA),cos(dphi)))
#
#
#    x <- c(cos(phiA[i])*cos(thetaA[i]),
#           sin(phiA[i]),
#           cos(phiA[i])*sin(thetaA[i]))
#    y <- R %*% x
#    print(R)
#    print(cbind(x,y))
#    phi[i] <- acos(y[3])
#    theta[i] <- acos(y[1]/sin(phi[i]))
#   }
  result <- list(phi=180*phi/pi,theta=180*theta/pi,)
  result
}
