# Computes Empirical Orthogonal Functions (EOFs)
#
# R.E. Benestad, met.no, Oslo, Norway 7.10.2002
# rasmus.benestad@met.no
#
#------------------------------------------------------------------------

EOF<-function(fields,l.wght=TRUE,lc180e=FALSE,direc="data/",
              lon=NULL,lat=NULL,l.stndrd=TRUE,las=1,
              mon=NULL,plot=TRUE,neofs=20,l.rm.ac=TRUE,lsave=TRUE,
              LINPACK=TRUE) {

#=========================================================================
library(ts)

if ((class(fields)[2]!="monthly.field.object") &
    (class(fields)[2]!="daily.field.object") &
    (class(fields)[1]!="field")) {
      print("class(fields) gives:")
      print(class(fields))
      stop("Need a 'field.object'")
    }

dir.0<-getwd()
if (!file.exists(direc)){
  print(paste("The directory",direc,"does not exists.. Creates it.."))
  dir.create(direc)
}

tim.unit <- attr(fields$tim,"unit")
tim.torg <- attr(fields$tim,"time_origin")

dims <- dim(fields$dat) 
if (length(dims)==3) dim(fields$dat) <- c(dims[1],dims[2]*dims[3])

# For naming the files containing the results

if (is.null(lon)) lon <- fields$lon
if (is.null(lat)) lat <- fields$lat
if (min(lon) < 0) deg.lon1.c<-"W" else deg.lon1.c<-"E"
if (max(lon) < 0) deg.lon2.c<-"W" else deg.lon2.c<-"E"
if (min(lat) < 0) deg.lat1.c<-"S" else deg.lat1.c<-"N"
if (max(lat) < 0) deg.lat2.c<-"S" else deg.lat2.c<-"N"
region<-paste(as.character(abs(round(min(lon)))),deg.lon1.c,
              as.character(abs(round(max(lon)))),deg.lon2.c,"-",
              as.character(abs(round(min(lat)))),deg.lat1.c,
              as.character(abs(round(max(lat)))),deg.lat2.c,sep="")
months<-c("Jan","Feb","Mar","Apr","May","Jun",
          "Jul","Aug","Sep","Oct","Nov","Dec")
season<-cbind(c(12,1,2),c(3,4,5),c(6,7,8),c(9,10,11))
season.c<-c("","DJF","MAM","JJA","SON")

id <- row.names(table(fields$id.x))
mm <- fields$mm
yy <- fields$yy
dd <- fields$dd
dat <- fields$dat
id.t <- fields$id.t
tim <- fields$tim
dims <- dim(dat)
nt <- dims[1]
np <- dims[2]
if (class(fields)[2]=="monthly.field.object") {
  if (is.null(mon))  {
      if (min(mm)==max(mm)) c.mon <- months[mm[1]]  else
               c.mon<-paste(months[min(mm)],"-",months[max(mm)],sep="")
      i.mm <- is.finite(mm)      
  } else {
      c.mon<-months[mon]
      i.mm <- is.element(mm,mon)
      dat <- dat[i.mm,]
      yy <- yy[i.mm]
      mm <- mm[i.mm]
      dd <- dd[i.mm]
      id.t <- id.t[i.mm]
      tim <- tim[i.mm]
    }
} else if (class(fields)[2]=="daily.field.object") {
  ac.mod<-matrix(rep(NA,nt*6),nt,6)
  ac.mod[,1]<-cos(2*pi*fields$tim/365.25)
  ac.mod[,2]<-sin(2*pi*fields$tim/365.25)
  ac.mod[,3]<-cos(4*pi*fields$tim/365.25)
  ac.mod[,4]<-sin(4*pi*fields$tim/365.25)
  ac.mod[,5]<-cos(6*pi*fields$tim/365.25)
  ac.mod[,6]<-sin(6*pi*fields$tim/365.25)
  if (l.rm.ac) {
    for (ip in seq(1,np,by=1)) {
      ac.fit<-lm(dat[,ip] ~ ac.mod)
      dat[!is.na(dat[,ip]),ip]<-ac.fit$residual
    }
  }
  if (is.null(mon))  {
    if (min(mm)==max(mm)) {
      c.mon <- months[mm[1]]
    } else { 
      c.mon<-paste(months[min(mm)],"-",months[max(mm)],sep="")
    }
    i.mm <- is.finite(mm) 
  } else {
    mon <- mod(mon-1,4)+1
    c.mon<-season.c[mon+1]
    mon <- season[,mon]
    i.mm <- is.element(mm,mon)
    dat <- dat[i.mm,]
    yy <- yy[i.mm]
    mm <- mm[i.mm]
    dd <- dd[i.mm]
    id.t <- id.t[i.mm]
    tim <- fields$tim[i.mm]
#    print(season[,mon])
  }
#  print("Months:")
#  print(table(mm))
#  print(paste("Season:",mon))
#  print(season)
} else if (substr(attr(fields$tim,"unit"),1,5)=="month") {
      c.mon<-months[as.numeric(row.names(table(mm)))]
      i.mm <- is.finite(mm) 
}

preds.names <- row.names(table(lower.case(fields$id.t)))
for (i.pred in 1:length(preds.names)) {
  if (preds.names[i.pred]=="t2m") preds.names[i.pred] <- tem
  if (preds.names[i.pred]=="psl") preds.names[i.pred] <- slp
}
preds.names <- row.names(table(preds.names))

preds.id <- ""
for (i.pred in 1:length(preds.names)) {
  eos <- nchar(preds.names[i.pred])

  if (instring("_",preds.names[i.pred])> 0) {
    eos <- instring("_",preds.names[i.pred])-1
    if (length(eos) > 1) eos <- eos[2]
  } else if (instring("-",preds.names[i.pred])> 0) {
    eos <- instring("-",preds.names[i.pred])-1
    if (length(eos) > 1) eos <- eos[2]
  } 
  preds.id  <- paste(preds.id,substr(preds.names[i.pred],1,eos),
                     "+",sep="")
}
print(preds.id)
preds.id <- substr(preds.id,1,nchar(preds.id)-1)
fname<-paste(direc,"eof_", preds.id,"_",fields$v.name,"_",region,"_",
       c.mon,'_',substr(attr(fields$tim,"unit"),1,3),".Rdata",sep="")
print(paste("File name:",fname,
            "sum(i.mm)=",sum(i.mm)))

#-------------------------------------------------------------------------

id <- row.names(table(fields$id.x))
size <- matrix(rep(0,3*fields$n.fld),3,fields$n.fld)
stdv <- rep(0,fields$n.fld)

#print(table(fields$id.lon))
#print(table(fields$id.lat))
#print(id)

for (i in 1:fields$n.fld) {
  ii <- fields$id.x == id[i]
#  print(paste("id.x: ",sum(ii)))
  id.lon <- fields$id.lon
  id.lat <- fields$id.lat
  i.lon <- id.lon == id[i]
  i.lat <- id.lat == id[i]
  lon.x <- fields$lon[i.lon]
  lat.x <- fields$lat[i.lat]
  id.lon <- id.lon[i.lon]
  id.lat <- id.lat[i.lat]
  nx <- length(lon.x)
  ny <- length(lat.x)
  dat.x <- dat[,ii]
  ix <- ((lon.x >= min(lon)) & (lon.x <= max(lon)))
  iy <- ((lat.x >= min(lat)) & (lat.x <= max(lat)))
  lon.x <- lon.x[ix]
  lat.x <- lat.x[iy]
  id.lon  <- id.lon[ix]
  id.lat  <- id.lat[iy]
#  print(dim(dat.x))
#  print(c(length(yy),ny,nx,sum(ix),sum(iy)))

  dim(dat.x) <- c(length(yy),ny,nx)
  dat.x <- dat.x[,iy,ix]

#print("Stdv[i]")  
  ny <- length(lat.x)
  nx <- length(lon.x)
  nt <- length(yy)
  stdv[i] <- sd(dat.x,na.rm=TRUE)

print("Remove mean values at each grid point")
  for (j.y in 1:ny) {
    for (i.x in 1:nx) {
      dat.x[,j.y,i.x] <- dat.x[,j.y,i.x] - mean(dat.x[,j.y,i.x],na.rm=TRUE)
    }
  }

  print("Add geographical weighting")
  if (l.wght) {
    print(paste("Weighting according to area. Field",i))
    Wght <-matrix(nrow=ny,ncol=nx)
    for (j in 1:nx)  Wght[,j]<-sqrt(abs(cos(pi*lat.x/180)))
    Wght[Wght < 0.01]<-NA     
    for (it in 1:nt) dat.x[it,,] <- dat.x[it,,]*Wght/stdv[i]
    print(paste("Wght.",i,"<-Wght",sep=""))
    eval(parse(text=paste("Wght.",i,"<-Wght",sep="")))
  }
  
  # reshape 3-D matrices to 2-D matrices

#print("Reshape 3-D matrices to 2-D matrices")
#print(dim(dat.x))
#print(c(nt,ny,nx,ny*nx))
  if (i == 1) {
    dat.d2 <- dat.x
    dim(dat.d2) <- c(nt,ny*nx)
  } else {
    dim(dat.x) <- c(nt,ny*nx)
    dat.d2 <- cbind(dat.d2,dat.x)
  }
  size[,i] <- c(nt,ny,nx)
  if (i==1) {
    id.x <- id[i]
    lons <- lon.x
    lats <- lat.x
    id.lons <- id.lon
    id.lats <- id.lat
#    print(paste("id.lons: ",length(id.lons),length(id.lon)))
#    print(paste("id.lats: ",length(id.lats),length(id.lat)))
  } else {
#    print(paste("Appending lons & lans, i=",i))
    id.x <- c(id.x,id[i])
    lons <- c(lons,lon.x)
    lats <- c(lats,lat.x)
    id.lons <- c(id.lons,id.lon)
    id.lats <- c(id.lats,id.lat)
#    print(paste("id.lons: ",length(id.lons),length(id.lon)))
#    print(paste("id.lats: ",length(id.lats),length(id.lat)))
  }
}

if (sum(is.na(dat.d2))>0) print(paste(sum(is.na(dat.d2)),
                            ' missing values of ',nx*ny*nt))
aver <- 0
print(paste("Find max autocorr in",ny*nx,"grid boxes."))
for (i in 1:(ny*nx)) {
  vec <- as.vector(dat.d2[,i])
  i.bad <- is.na(vec)
  if (sum(i.bad) == 0) {
    ar1 <- acf(vec[],plot=FALSE)
    aver <- max(c(aver,ar1$acf[2,1,1]))
  }
}
print(paste("aver=",aver))
n.eff <- round(nt * (1.0-aver)/(1.0+aver))  
print(paste("n.eff=",n.eff))

# Apply the PCA:       
print(paste("Singular Value Decomposition: ",sum(is.na(dat.d2)),
            ' NA-values -> set to zero of ',length(dat.d2)))
dat.d2[!is.finite(dat.d2)]<-0
print(range(dat.d2))
print(dim(dat.d2))
print(stdv)
if (LINPACK) pca<-svd(t(dat.d2)) else 
             pca<-La.svd(t(dat.d2))
PC<-pca$v[,1:neofs]
EOF<-t(pca$u[,1:neofs])
W<-pca$d[1:neofs]
tot.var <- sum(pca$d^2)
Var.eof<-100*pca$d[1:neofs]^2/tot.var

dW <- W*sqrt(2.0/n.eff)

# 2D->3D transform, invert weighting

#print("2D->3D transform")
i.last <- 0
x.last <- 0
y.last <- 0
for (i in 1:fields$n.fld) {
  i.fld <- seq(i.last+1,i.last+size[2,i]*size[3,i],by=1)
  i.last <- max(i.fld)
#  print(i.last)
  EOF.1 <- EOF[,i.fld]
  dim(EOF.1)<-c(neofs,size[2,i],size[3,i])
#  print(l.wght)
  if (l.wght) for (ieof in 1:neofs) EOF.1[ieof,,]<-
             EOF.1[ieof,,]*stdv[i]/eval(parse(text=paste("Wght.",i,sep="")))
#  print('eof.patt<-t(EOF.1[1,,])')
  eof.patt<-t(EOF.1[1,,])
  EOF[,i.fld] <- EOF.1
#  print('lonx,latx,...')
  lon.x <- lons[id.lons==id[i]]
  lat.x <- lats[id.lats==id[i]]
#  print('plot settings..')
  my.col <- rgb(c(seq(0,1,length=20),rep(1,21)),
                c(abs(sin((0:40)*pi/40))),
                c(c(rep(1,21),seq(1,0,length=20))))
  z.levs <- seq(-max(abs(as.vector(eof.patt)),na.rm=TRUE),
                max(abs(as.vector(eof.patt)),na.rm=TRUE),length=41)
#----------------------------------------------------------------
       
  if (plot) {
    filled.contour(lon.x,lat.x,eof.patt,levels=z.levs,
       main=paste("1st EOF for",id[i]),col=my.col,
       sub=paste(fields$f.name," (",c.mon,")"),
       xlab="Longitude",ylab="Latitude")
# From filled.contour in base
    mar.orig <- (par.orig <- par(c("mar","las","mfrow")))$mar
    on.exit(par(par.orig))

    w <- (3 + mar.orig[2]) * par('csi') * 2.54
    layout(matrix(c(2, 1), nc=2), widths=c(1, lcm(w)))
    
    par(las = las)
    mar <- mar.orig
    mar[4] <- 1
    par(mar=mar)
    contour(lon.x,lat.x,eof.patt,nlevels=7,add=TRUE,lwd=1,col="black")
    addland()
    grid()
    dev.copy2eps(file=paste("eof_",i,".eps",sep=""))
  }
}

attr(tim,"unit") <- tim.unit
attr(tim,"time_origin") <- tim.torg

#print("Construct list object")
eof<-list(EOF=EOF,W=W,PC=PC,id=preds.id,n.fld=fields$n.fld,tot.var=tot.var,
          id.t=id.t,id.x=fields$id.x,size=size,dW=dW,mon=mon,l.wght=l.wght,
          id.lon=id.lons,id.lat=id.lats,region=region,tim=tim,
          lon=lons,lat=lats,var.eof=Var.eof,yy=yy,mm=mm,dd=dd,
          v.name=fields$v.name,c.mon=c.mon,f.name=fname,
          attributes=fields$attributes)
class(eof) <- c("eof",class(fields))
save(file='data/ceof.Rdata',eof,ascii=FALSE)
if (lsave) save(file=fname,eof,ascii=FALSE) 

invisible(eof)
}
