# Optimal Interpolation procedure:
#
# After Reynolds and Smith (1994), J. Clim., June, pp 929--948
#
# R.E. Benestad

delta <- function(i,j) {
  if (i==j) delta <- 1 else delta <- 0
}

optint <- function(lon,lat,obs,
                   lon.grd,lat.grd,fguess,
                   eps,lambda=50,M=NULL,piipij=NULL,w=NULL,
                   tim=NULL,date=NULL) {

#  print('OPTINT: Optical Interpolation')
  a <- 6370  # Earth's radius
  np <- length(obs)
  nx <- length(lon.grd)
  ny <- length(lat.grd)
  dx <- 0.5*mean(diff(lon.grd))
  dy <- 0.5*mean(diff(lat.grd))
  obs.grd <- matrix(rep(NA,nx*ny),nx,ny)
  for (ip in 1:np) {
    ix <- (lon.grd >= lon[ip] - dx) &  (lon.grd <= lon[ip] + dx)
    iy <- (lat.grd >= lat[ip] - dy) &  (lat.grd <= lat[ip] + dy)
#    print(c(ip,lon.grd[ix],lat.grd[iy],obs[ip]))
    obs.grd[ix,iy] <- obs[ip]
  }

# Transform the coordinates from degrees to radians:
 
  if (is.null(piipij)) {  
    lon <- lon * pi/180
    lat <- lat * pi/180
    lons <- lon.grd
    lats <- lat.grd
    lon.grd <- lon.grd * pi/180
    lat.grd <- lat.grd * pi/180
#    print(paste('optint: Estimate the first-guess correlation-error',
#                '<pi pi>, dims=',ny*nx,'x',ny*nx,'=',ny^2*nx^2))
    piipij <- matrix(rep(0,ny^2*nx^2),ny*nx,ny*nx)
    for (j0 in 1:ny) {
      for (i0 in 1:nx) {
        ii <- (j0-1)*nx + i0
        for (j in 1:ny) {
          for (i in 1:nx) {
            jj <- (j-1)*nx + i
            piipij[ii,jj] <- exp( -( a^2*
              (lon.grd[i]*cos(lat.grd[j]) - lon.grd[i0]*cos(lat.grd[j0]))^2 +
              (lat.grd[j] - lat.grd[j0])^2 )/lambda^2 )
          }
        }
      }
    }
  }

#  x11()
#  image(piipij,main="piipij")
  
  if (is.null(M)) {
    dim(eps) <- c(nx*ny,1)
    M <- matrix(rep(0,ny^2*nx^2),ny*nx,ny*nx)
    for (j in 1:(ny*nx)) {
      for (i in 1:(ny*nx)) {
        M[i,j] <- piipij[i,j] + eps[i]*eps[j]*delta(i,j)
      }
    }
  }

#  x11()  
#  image(M,main="M")
#  x11()
#  image(M-piipij,main="M-piipij")
  
# Incomplete lower and upper triangular matrix decomposition:
# Solves the equation \sum_j M_ij w_ik = <pi_j pi_k> 

#  print("Solve 'sum_j{M_ij w_ik} = <pi_j pi_k>'")
  if (is.null(w)) {
    w <- matrix(rep(0,ny^2*nx^2),ny*nx,ny*nx)
    for (i in 1:ny*nx) {
      w[i,] <- qr.solve(M,piipij[,i]);
    }
  }
  
# Estimating the analysis increments: eq 1

  dim(obs.grd) <- c(nx*ny,1)
  dim(fguess) <- c(nx*ny,1)
  good <- !is.na(obs.grd)
  
  q <-  obs.grd - fguess
  r <- rep(0,ny*nx)
#  print("Solve 'r = w_i q'")

  for (i in 1:(ny*nx)) {
    r[i] <- w[i,good]*q[good]
  }
#  print(summary(as.vector(w)))
#  print(summary(q))
#  print(summary(r))
  dim(r) <- c(nx,ny)
  dim(fguess) <- c(nx,ny)

  results <- list(lon=lons,lat=lats,map=fguess + r,
                  tim=tim,date=date,M=M,piipij=piipij,w=w)
  class(results) <- "map"
#  attr(results) <- attr(obs)
  attr(results,"descr") <- "Optimal interpolation"
  invisible(results)
}
