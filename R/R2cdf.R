# Wrties an R/climpaxct-object as a CDF file that can
# be converted to netCDF through ncgen -b -o netcdf.file cdf.file

r2cdf <- function(filename,x,missing=-999.99,cleanup=TRUE,
                  ofs=NULL,scal=NULL) {

  if (class(x)=="field") {
    nt <- length(x$tim)
    ny <- length(x$lat)
    nx <- length(x$lon)
    neof <- NULL
  } else if (class(x)=="map") {
    nt <- 1
    ny <- length(x$lat)
    nx <- length(x$lon)
    neof <- NULL
  } else if (class(x)=="eof") {
    nt <- length(x$tim)
    ny <- length(x$lat)
    nx <- length(x$lon)
    neof <- length(x$lam)
  } else {
    print("x has none of the classes: field, map, or EOF")
    return()
  }

  if ((is.null(ofs)) & (min(x$dat,na.rm=TRUE) > 0)) {
    ofs <- 10^(round(log(min(x$dat,na.rm=TRUE))/log(10)))
  }
  if (is.null(scal)) {
    max.dev <- range(x$dat[is.finite(x$dat)])
    scal <- 10^(-round(log(32000/(max.dev[2]-max.dev[1]))/log(10)))
  }

  x$dat[!is.finite(x$dat)] <- missing
  cdf <- file(paste(filename,".cdf",sep=""),"w")
  cat("netcdf DNMI_slp {","dimensions:",file=cdf,sep = "\n")
  cat(paste("        Lon =",nx,";"),file=cdf,sep = "\n")
  cat(paste("        Lat =",ny,";"),file=cdf,sep = "\n")
  cat(paste("        Time =",nt,";"),file=cdf,sep = "\n")
  if (!is.null(neof)) cat(paste("eof =",neof,";"),file=cdf,sep = "\n")

  cat("variables:",file=cdf,sep = "\n")
  cat("        float Lon(Lon) ;",file=cdf,sep = "\n")
  cat('                Lon:units = "degrees_east" ;',file=cdf,sep = "\n")
  cat('                Lon:modulo = " " ;',file=cdf,sep = "\n")
  cat('                Lon:long_name = "longitude" ;',file=cdf,sep = "\n")
  cat("        float Lat(Lat) ;",file=cdf,sep = "\n")
  cat('                Lat:units = "degrees_north" ;',file=cdf,sep = "\n")
  cat('                Lat:long_name = "latitude" ;',file=cdf,sep = "\n")
  cat("        float Time(Time) ;",file=cdf,sep = "\n")
  cat(paste('                Time:units = "',x$attributes$time.unit,'" ;',sep=""),
      file=cdf,sep = "\n")
  cat(paste('                Time:time_origin = "',x$attributes$time.origin,
            '" ;',sep=""),file=cdf,sep = "\n")
  if (!is.null(neof)) {
    cat("        short eof(eof) ;",file=cdf,sep = "\n")
    cat('                eof:units = " " ;',file=cdf,sep = "\n")
    cat('                eof:long_name = "mode" ;',file=cdf,sep = "\n")
  }

  if (is.null(neof)) {
    cat(paste("        short ",x$v.name,"(Time, Lat, Lon) ;",sep=""),
        file=cdf,sep = "\n")
    cat(paste('              ',x$v.name,':units = "',x$attributes$unit,'" ;',sep=""),
        file=cdf,sep = "\n")
    cat(paste('              ',x$v.name,':longname = "',
              x$attributes$longname,'" ;',sep=""),
        file=cdf,sep = "\n")
    cat(paste('              ',x$v.name,':missing_value = ',missing,' ;',sep=""),
        file=cdf,sep = "\n")
    cat(paste('              ',x$v.name,':add_offset = ',ofs,' ;',sep=""),
        file=cdf,sep = "\n")
        cat(paste('              ',x$v.name,':scale_factor = ',scal,' ;',sep=""),
        file=cdf,sep = "\n")
  } else {
    cat("        float EOF(neof,Lat,Lon) ;",file=cdf,sep = "\n")
    cat('                EOF:units = " " ;',file=cdf,sep = "\n")
    cat('                EOF:long_name = "Empirical Orthogonal Functions" ;',
        file=cdf,sep = "\n")
    cat("        float PC(Time,neof) ;",file=cdf,sep = "\n")
    cat('                PC:units = " " ;',file=cdf,sep = "\n")
    cat('                PC:long_name = "Principal Components" ;',
        file=cdf,sep = "\n")
    cat("        float lambda(neof) ;",file=cdf,sep = "\n")
    cat('                lambda:units = " " ;',file=cdf,sep = "\n")
    cat('                lambda:long_name = "Eigenvalues" ;',
        file=cdf,sep = "\n")
  }
    
  cat(" ",file=cdf,sep = "\n")
  cat("// global attributes:",file=cdf,sep = "\n")
  cat('                   :history = "Saved by r2cdf.R (clim.pact)" ;',
      file=cdf,sep = "\n")
  cat('                   :URL = "http://cran.r-project.org/" ;',
      file=cdf,sep = "\n")
      
  cat("data:",file=cdf,sep = "\n")
  cat(" ",file=cdf,sep = "\n")

  cat(" Lon = ",file=cdf,sep = " ")
  cat(as.character(x$lon[1:(nx-1)]),file=cdf,sep = ", ")
  cat(as.character(x$lon[nx]),file=cdf,sep = " ")
  cat(";",file=cdf,sep = "\n")
  cat(" ",file=cdf,sep = "\n")
  
  cat(" Lat = ",file=cdf,sep = " ")
  cat(as.character(x$lat[1:(ny-1)]),file=cdf,sep = ", ")
  cat(as.character(x$lat[ny]),file=cdf,sep = " ")
  cat(";",file=cdf,sep = "\n") 
  cat(" ",file=cdf,sep = "\n")
  
  cat(" Time = ",file=cdf,sep = " ")
  cat(as.character(x$tim[1:(nt-1)]),file=cdf,sep = ", ")
  cat(as.character(x$tim[nt]),file=cdf,sep = " ")
  cat(";",file=cdf,sep = "\n") 
  cat(" ",file=cdf,sep = "\n")
  
  if (!is.null(neof)) {
    cat(" eof = ",file=cdf,sep = " ")
    cat(as.character(1:(neof-1)),file=cdf,sep = ", ")
    cat(as.character(neof),file=cdf,sep = " ")
    cat(";",file=cdf,sep = "\n")
    cat(" ",file=cdf,sep = "\n")
  }

  if (is.null(neof)) {
    cat(paste(" ",x$v.name,sep=""),file=cdf,sep = " ")
    cat("= ",file=cdf,sep = "\n")
    for (it in 1:nt) {
      for (j in 1:ny) {
        cat(as.character(round((x$dat[it,j,]-ofs)/scal)),file=cdf,sep = ", ")
        if ((it < nt) | (j < ny)) cat(", ",file=cdf,sep = "\n") else
                                  cat("; ",file=cdf,sep = "\n")
      }
    }
  } else {
    cat(" EOF = ",file=cdf,sep = "\n")
    for (it in 1:neof) {
      for (j in 1:ny) {
        cat(as.character(round(x$EOF[it,j,])),file=cdf,sep = " ")
      }
    }
    cat(" PC = ",file=cdf,sep = "\n")
    for (it in 1:neof) cat(as.character(round(x$PC[,it]),4),file=cdf,sep = " ")
    cat(" Lambda = ",file=cdf,sep = "\n")
    cat(as.character(round(x$lambda,4)),file=cdf,sep = " ")
  }
  cat("}",file=cdf,sep = "\n")
  close(cdf)
  
  system(paste("ncgen -b  -o ",filename,".nc ",filename,".cdf",sep=""),intern=T)

  if (cleanup) system(paste("rm -f ",filename,sep=""),intern=T)
}
