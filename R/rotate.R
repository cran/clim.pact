rotate <- function(lons,lats,lon.0=NULL,lat.0=NULL,method="Cayley-Klein",test=TRUE) {

# lons.0 and lats.0 is the new central coordinate...

  if (is.null(lon.0)) lon.0 <- mean(lons)
  if (is.null(lat.0)) lat.0 <- mean(lats)

  theta0 <- pi*lon.0/180
  phi0 <- pi*lat.0/180
  r0 <- c(cos(phi0)*cos(theta0),
          sin(phi0),
          cos(phi0)*sin(theta0))

  theta <- rep(NA,length(lons)); phi <- theta

  if (method=="old") {
  
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
  } else {

# URL: http://www.uwgb.edu/dutchs/mathalgo/sphere0.htm
# http://www.iop.org/EJ/article/0305-4470/32/26/101/ja32026l1.html
# Let the coordinates of the rotation axis be (c1,c2,c3) and the rotation angle be 'a'. 
# Any point (x,y,z) will be rotated into new coordinates (x',y',z') by the following matrix:
#
#     |1 0 0|              |c1c1 c1c2 c1c3|        |  0 -c3  c2|
#cos a|0 1 0| + (1 - cos a)|c2c1 c2c2 c2c3| + sin a| c3   0 -c1| 
#     |0 0 1|              |c3c1 c3c2 c3c3|        |-c2  c1   0|
#
#In the middle term the coefficients are written c1c1, c2c2 and c3c3 partly to avoid exponents but 
#also to show the pattern of the terms more clearly.
#
# k = c1 sin a/2, m = c2 sin a/2, n = c3 sin a/2, p = cos a/2. Note that:
#
#                 2    2    2    2
#                k  + m  + n  + p  = 1
#
#The rotation matrix is:
#
#          | 2  2  2  2                         | 
#          |k -m -n +p    2(km-np)    2(nk+mp)  |
#          |                                    |
#          |              2  2  2  2            | 
#          | 2(km-np)    m -m -m +p   2(mn-kp)  |
#          |                                    |
#          |                          2  2  2  2| 
#          | 2(nk-mp)     2(mn+kp)   m -m -m +p |
#
# c1=cos(lat)cos(long)
# c2=cos(lat)sin(long)
# c3=sin(lat)
# x' = x cos a + (1 - cos a)(c1c1x + c1c2y + c1c3z) + (c2z - c3y)sin a
# y' = y cos a + (1 - cos a)(c2c1x + c2c2y + c2c3x) + (c3x - c1z)sin a
# z' = z cos a + (1 - cos a)(c3c1x + c3c2y + c3c3z) + (c1y - c2x)sin a

     good <- is.finite(lons) & is.finite(lats)
     lons <- lons[good]; lats <- lats[good]
     a1 <- theta0
     a2 <- phi0

#print(c(a1,a2))

       thetaA <- pi*lons/180
       phiA <- pi*lats/180
       r <- rbind(cos(phiA)*cos(thetaA),
                  cos(phiA)*sin(thetaA),
                  sin(phiA))       

# longitudal rotation:

       c1 <- 0; c2 <- 0; c3 <- 1
       k <- c1 * sin(a1/2)
       m <- c2 * sin(a1/2)
       n <- c3 * sin(a1/2)
       p <- cos(a1/2)

       rot.mat <- matrix(c(k^2-m^2-n^2+p^2,2*(k*m-n*p),2*(n*k+m*p),
                           2*(k*m-n*p),p^2-m^2,2*(m*n-k*p),
                           2*(n*k-m*p),2*(m*n+k*p),p^2-m^2),
                           nrow=3,ncol=3,byrow = FALSE)


#print(round(r,2))
#print(rot.mat)
       r <- rot.mat %*% r
#print(round(r,2))
       R <- colSums(r)

# latitudal rotation:

       if (a2 > 0) {
       c1 <- -R[2]/sqrt(R[1]^2 + R[2]^2); c2 <- R[1]/sqrt(R[1]^2 + R[2]^2); c3 <- 0
       k <- c1 * sin(a2/2)
       m <- c2 * sin(a2/2)
       n <- c3*sin(a2/2)
       p <- cos(a2/2)

       rot.mat <- matrix(c(k^2-m^2-n^2+p^2,2*(k*m-n*p),2*(n*k+m*p),
                           2*(k*m-n*p),p^2-m^2,2*(m*n-k*p),
                           2*(n*k-m*p),2*(m*n+k*p),p^2-m^2),
                           nrow=3,ncol=3,byrow = TRUE)

       r <- rot.mat %*% r
       }

       phi <-   asin( r[3,] )
       theta <- acos(r[1,]/cos(phi))
       phi <- 180*phi/pi; theta <- 180*theta/pi
       theta[theta > 180] <- theta - 360
     
   if (test) {
      print("test")
      plot(lons,lats,pch=".")
      points(theta,phi,pch=".",col="grey")
      print(summary(theta))
   }
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
