
load("eof_DNMI_slp.Rdata")

  xxx <-t(eof$EOF) %*% diag(eof$W) %*% t(eof$PC); xxx <- t(xxx)
  dim(xxx) <- eof$size
  climxxx <- t(matrix(eof$clim,eof$size[3],eof$size[2]))
  for (i in 1:eof$size[1]) xxx[i,,] <- xxx[i,,] + climxxx
  DNMI.slp  <- list(dat=xxx,lon=eof$lon,lat=eof$lat,tim=eof$tim,lev=NULL,
                       v.name=eof$v.name,id.x=eof$id.x,id.t=eof$id.t,
                       yy=eof$yy,mm=eof$mm,dd=eof$dd,n.fld=eof$n.fld,
                       id.lon=eof$id.lon,id.lat=eof$id.lat,
                       attributes=eof$attributes,filename="DNMI.slp")
  class(DNMI.slp) <- class(eof)[-1]
rm(xxx,climxxx,i)
