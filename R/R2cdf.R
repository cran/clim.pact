# Writes an R/climpaxct-object as a CDF file that can
# be converted to netCDF through ncgen -b -o netcdf.file cdf.file
# R.E. Benestad

r2cdf <- function(filename,x,missing=-999.99,cleanup=TRUE,
                  ofs=NULL,scal=NULL) {

  if (class(x)[1]=="field") {
    nt <- length(x$tim)
    ny <- length(x$lat)
    nx <- length(x$lon)
    neof <- NULL
  } else if (class(x)[1]=="map") {
    nt <- 1
    ny <- length(x$lat)
    nx <- length(x$lon)
    x$dat <- t(x$map); dim(x$dat) <- c(nt,ny,nx)
    if (is.null(x$v.name)) x$v.name <- "map"
    if (is.null(x$attributes)) {
      x$attributes <- list(time.unit="unknown",
                           time.origin="unknown",
                           unit="unknown",
                           longname="unknown")
    }
    neof <- NULL
    x$map[!is.finite(x$map)] <- missing
    if (is.null(x$tim)) x$tim <- 0
  } else if (class(x)[1]=="eof") {
    nt <- length(x$tim)
    ny <- length(x$lat)
    nx <- length(x$lon)
    x$dat <- x$EOF
    id <- row.names(table(x$id.x))
    neof <- length(x$W)
  } else {
    print("x has none of the classes: field, map, or EOF")
    return()
  }
  
  if ((is.null(ofs)) & (min(x$dat,na.rm=TRUE) > 0) & (class(x)[1]=="field")) {
    ofs <- 10^(round(log(min(x$dat,na.rm=TRUE))/log(10)))
  } else ofs <- 0
  if ((is.null(scal)) & (class(x)[1]=="field")) {
    max.dev <- range(x$dat[is.finite(x$dat)])
    scal <- 10^(-round(log(32000/(max.dev[2]-max.dev[1]))/log(10)))
  } else scal <- 1

  x$dat[!is.finite(x$dat)] <- missing
  cdf <- file(paste(filename,".cdf",sep=""),"w")
  cat("netcdf DNMI_slp {","dimensions:",file=cdf,sep = "\n")
  if (class(x)[1]!="eof") {
    cat(paste("        Lon =",nx,";"),file=cdf,sep = "\n")
    cat(paste("        Lat =",ny,";"),file=cdf,sep = "\n")
  } else {
    for (i in 1:x$n.fld) {
      i.lon <- x$id.lon == id[i]
      i.lat <- x$id.lat == id[i]
      lon.x <- x$lon[i.lon]
      lat.x <- x$lat[i.lat]
      nx <- length(lon.x)
      ny <- length(lat.x)
      cat(paste("        Lon",i," = ",nx," ;",sep=""),file=cdf,sep = "\n")
      cat(paste("        Lat",i," = ",ny," ;",sep=""),file=cdf,sep = "\n")
    }
  }
  cat(paste("        Time =",nt,";"),file=cdf,sep = "\n")
  if (!is.null(neof)) {
    cat(paste("        eof =",neof,";"),file=cdf,sep = "\n")
  }

  cat("variables:",file=cdf,sep = "\n")
  if (class(x)[1]!="eof") {
    cat("        float Lon(Lon) ;",file=cdf,sep = "\n")
    cat('                Lon:units = "degrees_east" ;',file=cdf,sep = "\n")
    cat('                Lon:modulo = " " ;',file=cdf,sep = "\n")
    cat('                Lon:long_name = "longitude" ;',file=cdf,sep = "\n")
    cat("        float Lat(Lat) ;",file=cdf,sep = "\n")
    cat('                Lat:units = "degrees_north" ;',file=cdf,sep = "\n")
    cat('                Lat:long_name = "latitude" ;',file=cdf,sep = "\n")
  } else {
    for (i in 1:x$n.fld) {
      i.lon <- x$id.lon == id[i]
      i.lat <- x$id.lat == id[i]
      lon.x <- x$lon[i.lon]
      lat.x <- x$lat[i.lat]
      nx <- length(lon.x)
      ny <- length(lat.x)
      cat(paste("        float Lon",i,"(Lon",i,") ;",sep=""),
          file=cdf,sep = "\n")
      cat(paste('                Lon',i,':units = "degrees_east" ;',sep=""),
          file=cdf,sep = "\n")
      cat(paste('                Lon',i,':modulo = " " ;',sep=""),file=cdf,
          sep = "\n")
      cat(paste('                Lon',i,':long_name = "longitude" ;',sep=""),
          file=cdf,sep = "\n")
      cat(paste("        float Lat",i,"(Lat",i,") ;",sep=""),file=cdf,
          sep = "\n")
      cat(paste('                Lat',i,':units = "degrees_north" ;',sep=""),
          file=cdf,sep = "\n")
      cat(paste('                Lat',i,':long_name = "latitude" ;',sep=""),
          file=cdf,sep = "\n")
    }
  }
  
  cat("        float Time(Time) ;",file=cdf,sep = "\n")
  cat(paste('                Time:units = "',x$attributes$time.unit,'s since ',x$attributes$time.origin,'" ;',
          sep=""),file=cdf,sep = "\n")
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
    for (i in 1:x$n.fld) {
      cat(paste("        float EOF",i,"(eof,Lat",i,",Lon",i,") ;",sep=""),
          file=cdf,sep = "\n")
      cat(paste('                EOF',i,':units = " " ;',sep=""),
          file=cdf,sep = "\n")
      cat(paste('                EOF',i,
                ':long_name = "Empirical Orthogonal Functions" ;',sep=""),
        file=cdf,sep = "\n")
    }
    cat("        float PC(Time,eof) ;",file=cdf,sep = "\n")
    cat('                PC:units = " " ;',file=cdf,sep = "\n")
    cat('                PC:long_name = "Principal Components" ;',
        file=cdf,sep = "\n")
    cat("        float lambda(eof) ;",file=cdf,sep = "\n")
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

  if (class(x)[1]!="eof") {
    cat(" Lon = ",file=cdf,sep = " ")
    cat(as.character(round(x$lon,6)),file=cdf,sep = ", ")
    cat(";",file=cdf,sep = "\n")
    cat(" ",file=cdf,sep = "\n")
  
    cat(" Lat = ",file=cdf,sep = " ")
    cat(as.character(round(x$lat,6)),file=cdf,sep = ", ")
    cat(";",file=cdf,sep = "\n") 
    cat(" ",file=cdf,sep = "\n")
  } else {
    for (i in 1:x$n.fld) {
      i.lon <- x$id.lon == id[i]
      i.lat <- x$id.lat == id[i]
      lon.x <- x$lon[i.lon]
      lat.x <- x$lat[i.lat]
      nx <- length(lon.x)
      ny <- length(lat.x)
      cat(paste(" Lon",i," = ",sep=""),file=cdf,sep = " ")
      cat(as.character(round(lon.x,6)),file=cdf,sep = ", ")
      cat(";",file=cdf,sep = "\n")
      cat(" ",file=cdf,sep = "\n")
  
      cat(paste(" Lat",i," = ",sep=""),file=cdf,sep = " ")
      cat(as.character(round(lat.x,6)),file=cdf,sep = ", ")
      cat(";",file=cdf,sep = "\n") 
      cat(" ",file=cdf,sep = "\n")
    }
}
  cat(" Time = ",file=cdf,sep = " ")
  if (class(x)[1]!="map"){
#    cat("0",file=cdf,sep = ", ")
    cat(as.character(round(x$tim,1)),file=cdf,sep = ", ")
    cat(";",file=cdf,sep = "\n") 
    cat(" ",file=cdf,sep = "\n")
  } else {
    cat(as.character(x$tim),";",file=cdf,sep = "\n")
    cat(" ",file=cdf,sep = "\n")
  }
  if (!is.null(neof)) {
    cat(" eof = ",file=cdf,sep = " ")
    cat(as.character(1:neof),file=cdf,sep = ", ")
    cat(";",file=cdf,sep = "\n")
    cat(" ",file=cdf,sep = "\n")
  }

  if (class(x)[1]=="map"){
    cat(paste(" ",x$v.name,sep=""),file=cdf,sep = " ")
    cat("= ",file=cdf,sep = "\n")
    for (j in 1:ny) {
      cat(as.character(round((x$map[,j]-ofs)/scal)),file=cdf,sep = ", ")
      if (j < ny) cat(", ",file=cdf,sep = "\n") else
                  cat("; ",file=cdf,sep = "\n")
    }
  } else   if (class(x)[1]=="field"){
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
    i.last <- 0
    for (i in 1:x$n.fld) {
#      print(i)
#      print(x$size[,i])
      i.fld <- seq(i.last+1,i.last+x$size[2,i]*x$size[3,i],by=1)
      i.last <- max(i.fld)
      EOF.1 <- x$EOF[,i.fld]
#      print(dim(EOF.1))
#      print(c(neof,ny,nx))
      
      dim(EOF.1)<-c(neof,x$size[2,i],x$size[3,i])
      EOF.1[!is.finite(EOF.1)] <- missing
      i.lon <- x$id.lon == id[i]
      i.lat <- x$id.lat == id[i]
      lon.x <- x$lon[i.lon]
      lat.x <- x$lat[i.lat]
      nx <- length(lon.x)
      ny <- length(lat.x)
      
      cat(paste(" EOF",i," = ",sep=""),file=cdf,sep = "\n")
      for (it in 1:neof) {
        for (j in 1:ny) {
          cat(as.character(round(EOF.1[it,j,],3)),file=cdf,sep = ",  ")
          if ((it < neof) | (j < ny)) cat(", ",file=cdf,sep = "\n") else
                                      cat("; ",file=cdf,sep = "\n")
        }
      }
      cat(" ",file=cdf,sep = "\n")
      
      cat(" PC = ",file=cdf,sep = "\n")
      for (it in 1:neof) {
        cat(as.character(round(x$PC[,it],3)),file=cdf,sep = ",  ")
        if (it < neof) cat(", ",file=cdf,sep = "\n") else
                       cat("; ",file=cdf,sep = "\n")
      }
      cat(" ",file=cdf,sep = "\n")
      
      cat(" lambda = ",file=cdf,sep = "\n")
      cat(as.character(round(x$W,4)),file=cdf,sep = ",  ")
      cat("; ",file=cdf,sep = "\n")
    }
  }
  cat("}",file=cdf,sep = "\n")
  close(cdf)
  
  system(paste("ncgen -b  -o ",filename,".nc ",filename,".cdf",sep=""),intern=T)

  if (cleanup) system(paste("rm -f ",filename,sep=""),intern=T)
}
