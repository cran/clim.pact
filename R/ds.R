# Empirical downscaling using EOFs of monthly values from eof.R
# Predictand is a time series of monthly values from NACD or climate station.
#
# Reference: R.E. Benestad et al. (2002),
#            Empirically downscaled temperature scenarios for Svalbard,
#            doi.10.1006/asle.2002.005, September 18.
#
#            R.E. Benestad (2001),
#            A comparison between two empirical downscaling strategies,
#            Int. J. Climatology, 1645-1668, vol. 21, DOI 10.1002/joc.703
#
# R.E. Benestad, met.no, Oslo, Norway 16.04.2002
# rasmus.benestad@met.no
#------------------------------------------------------------------------

DS <- function(dat,preds,mon=NULL,direc="output/",cal.id=NULL,
               ldetrnd=TRUE,i.eofs=seq(1,8,by=1),ex.tag="",
               method="lm",plot=TRUE,leps=FALSE,param="t2m",
               plot.res=FALSE,plot.rate=FALSE,xtr.args="",
               swsm="step",predm="predict",lsave=TRUE,rmac=TRUE,
               silent=FALSE) {
library(ts)
library(ctest)
#library(chron)
#library(date)
library(xtable)

dir.0<-getwd()
if (!file.exists(direc)){
  if (!silent) print(paste("The directory",direc,"does not exists.. Creates it.."))
  dir.create(direc)
}
if (class(preds)[1]!="eof") {
  stop("The predictor must be an 'eof' object.")
}

if (class(dat)[1]!="station") {
  stop(paste("The predictand must be a 'monthly.station.record'",
             "object - Use  station.obj()"))
}

anm.weight <- FALSE
if (method=="anm.weight") {
  method <- "anm"
  anm.weight <- TRUE
}
if (method=="anm") {
  swsm <- "none"
#  predm <- "predictANM"
  ldetrnd <- FALSE
  rmac <- FALSE
}

if (class(preds)[2]=="daily.field.object") {
  good <- eval(parse(text=paste("is.finite(dat$",param,")",sep="")))
  eval(parse(text=paste("dat$",param," <- dat$",param,"[good]",sep="")))
  dat$mm <- dat$mm[good]; dat$yy <- dat$yy[good]
  dat$dd <- dat$dd[good]; dat$tim <- dat$tim[good]
} 

cmon<-c('Jan','Feb','Mar','Apr','May','Jun',
        'Jul','Aug','Sep','Oct','Nov','Dec')
season.c<-c("","DJF","MAM","JJA","SON")
season<-cbind(c(12,1,2),c(3,4,5),c(6,7,8),c(9,10,11))
lon <- preds$lon
lat <- preds$lat
if (min(lon) < 0) deg.lon1.c<-"W" else deg.lon1.c<-"E"
if (max(lon) < 0) deg.lon2.c<-"W" else deg.lon2.c<-"E"
if (min(lat) < 0) deg.lat1.c<-"S" else deg.lat1.c<-"N"
if (max(lat) < 0) deg.lat2.c<-"S" else deg.lat2.c<-"N"
region<-paste(as.character(abs(round(min(lon)))),deg.lon1.c,
              as.character(abs(round(max(lon)))),deg.lon2.c,"-",
              as.character(abs(round(min(lat)))),deg.lat1.c,
              as.character(abs(round(max(lat)))),deg.lat2.c,sep="")

month <-cmon[mon]
if ((class(preds)[2]=="daily.field.object") & !is.null(mon)) {
  mon <- mod(mon-1,4)+1
  month <- season.c[mon+1]
  mon <- season[mon]
}

if (!is.null(mon) & !is.null(preds$mon)) {    
  if (is.null(mon)) mon <- preds$mon
  if (!silent) print(paste("Extract",cmon[mon],"-> # of data points=",
              sum(is.element(preds$mm,mon))))
  if (  sum(is.element(preds$mm,mon))==0  ) {
    if (!silent) print(paste(">>> ",cmon[mon],
                             " is not found in the PCA product! <<<"))
    months <- row.names(table(preds$mm))
    if (!silent) print(paste(" Available months are:",cmon[months]))
    mon<-months
  }
} else {
  if (max(preds$mm) > min(preds$mm)) {
    month <- paste(cmon[min(preds$mm)],"-",cmon[max(preds$mm)],sep="")
  } else month  <- cmon[mean(preds$mm)]
  months <- row.names(table(preds$mm))
  mon<-months
}

preds.names <- row.names(table(preds$id.t))
preds.id <- ""
for (i.pred in 1:length(preds.names)) {
  eos <- nchar(preds.names[i.pred])
  if (instring("-",preds.names[i.pred])> 0) {
    eos <- instring("-",preds.names[i.pred])-1
  } else if (instring("_",preds.names[i.pred])> 0) {
    eos <- instring("_",preds.names[i.pred])-1
  }
  preds.id  <- paste(preds.id,substr(preds.names[i.pred],1,eos),
                     "+",sep="")
}
if (is.null(attr(preds$tim,"unit"))) {
  if (preds$dd[2]-preds$dd[1]==0) attr(preds$tim,"unit")<-"mon" else
                                 attr(preds$tim,"unit")<-"day"
#  if (!silent) print(attr(preds$tim,"unit"))
}

eos <- instring(" ",dat$location)[1]-1
if ((is.null(eos)) | (eos <= 0)) eos <- nchar(dat$location)
preds.id <- substr(preds.id,1,nchar(preds.id)-1)
fname<-paste(direc,"ds_",preds.id,"_",region,"_",
             substr(dat$location,1,eos),"_",dat$ele,"_",preds$c.mon,'_',
             substr(attr(preds$tim,"unit"),1,3),"_",method,
             ex.tag,".Rdata",sep="")

# Get the predictand

loc <- dat$location
if (class(dat)[2]=="monthly.station.record"){
  v.name  <- abbreviate(dat$obs.name)
} else if (class(dat)[2]=="daily.station.record") {
  v.name  <- param
}

if (v.name=="mT(2") v.name <- "T"
ny<-length(dat$yy)

if ( (ny < 20) | (sum(is.na(dat$yy))>0) ) {
  if (!silent) print("ds: WARNING: ... SENSING POSSIBLE PROBLEMS!...")
  if (!silent) print(paste("For predictor, you selected",preds$f.name))
  if (!silent) print(paste("for predictand, you selected",loc))
  if (!silent) print(paste("Number of valid data points from station",
              sum(!is.na(dat$val)),"- Length of data record=",ny))
  if (!silent) print("Years of station observations:")
  if (!silent) print(range(dat$yy))
  if (!silent) print("Years of predictor (PCs):")
  if (!silent) print(range(preds$yy))
  if (!silent) print("You may want to try another location or different dataset")
}

if (class(dat)[2]=="monthly.station.record") {
  yy.o<-sort(rep(dat$yy,12))
  mm.o<-rep(seq(1,12,by=1),ny)
  dd.o <- rep(15,length(yy.o))
  y.o<-t(dat$val)
  dim(y.o)<-c(12*ny,1)
  ds.unit <- dat$unit
} else if (class(dat)[2]=="daily.station.record") {
  yy.o <- dat$yy
  mm.o <- dat$mm
  dd.o <- dat$dd
  if (eval(parse(text=paste("is.null(dat$",param,")",sep="")))) {
    if (!silent) print(summary(dat))
    param <- readline("Select object field:")
    ds.unit <- readline("unit:")
  } else {
    if (!silent) print(paste("y.o<-dat$",param,sep=""))
    eval(parse(text=paste("y.o<-dat$",param,sep="")))
    ds.unit <- dat$unit[1]
  }
#  tim.o <- julian(mm.o,dd.o,yy.o,origin.=c(1,1,1970))
#  tim.o <- mdy.date(mm.o, dd.o, yy.o)
  tim.o <- julday(mm.o, dd.o, yy.o) - julday(1,1,1970)
  nt <- length(tim.o)
  ac.mod<-matrix(rep(NA,nt*6),nt,6)
  ac.mod[,1]<-cos(2*pi*tim.o/365.25); ac.mod[,2]<-sin(2*pi*tim.o/365.25)
  ac.mod[,3]<-cos(4*pi*tim.o/365.25); ac.mod[,4]<-sin(4*pi*tim.o/365.25)
  ac.mod[,5]<-cos(6*pi*tim.o/365.25); ac.mod[,6]<-sin(6*pi*tim.o/365.25)
  ac.obs <- data.frame(y=y.o, X=ac.mod)
  ac.fit<-lm(y ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6,data=ac.obs)
  if (rmac) y.o <- ac.fit$residual
}

y.o[y.o < -99] <-NA
mm.o <- mm.o[!is.na(y.o)]
yy.o <- yy.o[!is.na(y.o)]
dd.o <- dd.o[!is.na(y.o)]
y.o <- y.o[!is.na(y.o)]

# Give the correct file names.
n.fld <- preds$n.fld
if (n.fld>1) {
  reg <- paste(preds$v.name[1],preds$v.name[2],"_",
               as.character(-1*min(preds$lon)),"W-",
               as.character(max(preds$lon)),"E_",
               as.character(min(preds$lat)) ,"N-",
               as.character(max(preds$lat)),"N",sep="")
} else {
  reg <- paste(preds$v.name[1],"_",
               as.character(-1*min(preds$lon)),"W-",
               as.character(max(preds$lon)),"E_",
               as.character(min(preds$lat)) ,"N-",
               as.character(max(preds$lat)),"N",sep="")
}

y.o<-  y.o[is.element(mm.o,mon)]
yy.o<- yy.o[is.element(mm.o,mon)]
dd.o<- dd.o[is.element(mm.o,mon)]
mm.o<- mm.o[is.element(mm.o,mon)]

if (is.null(cal.id)) cal.id<- preds$id.t[1]
#print("Calibration predictors:")
#print(sum(preds$id.t==cal.id & !is.na(preds$PC[,1])))
#print(range(preds$yy[preds$id.t==cal.id & !is.na(preds$PC[,1])]))
X.cal<-  preds$PC[preds$id.t==cal.id & !is.na(preds$PC[,1]),]
yy.cal<- preds$yy[preds$id.t==cal.id & !is.na(preds$PC[,1])]
mm.cal<- preds$mm[preds$id.t==cal.id & !is.na(preds$PC[,1])]
dd.cal<- preds$dd[preds$id.t==cal.id & !is.na(preds$PC[,1])]

if (!silent) print("------------Match times---------------- ")
X.cal<-  X.cal[is.element(mm.cal,mon),]
yy.cal<- yy.cal[is.element(mm.cal,mon)]
dd.cal<- dd.cal[is.element(mm.cal,mon)]
mm.cal<- mm.cal[is.element(mm.cal,mon)]
#print(c(length(y.o),length(yy.o),length(X.cal[,1]),length(yy.cal)))

if (sum((preds$id.t!=cal.id) & !is.na(preds$PC[,1]))>0) {
  X.gcm<-preds$PC[preds$id.t!=cal.id  & !is.na(preds$PC[,1]),]
  yy.gcm<-preds$yy[preds$id.t!=cal.id & !is.na(preds$PC[,1])]
  mm.gcm<-preds$mm[preds$id.t!=cal.id & !is.na(preds$PC[,1])]
  dd.gcm<-preds$dd[preds$id.t!=cal.id & !is.na(preds$PC[,1])]
  
#print("The scenarios:")
  X.gcm<-X.gcm[is.element(mm.gcm,mon),]
  yy.gcm<-yy.gcm[is.element(mm.gcm,mon)]
  dd.gcm<-dd.gcm[is.element(mm.gcm,mon)]
  mm.gcm<-mm.gcm[is.element(mm.gcm,mon)]
} else {
  X.gcm<-X.cal
  yy.gcm <- yy.cal
  mm.gcm <- mm.cal
  dd.gcm <- dd.cal
}

# Find the common period:

if (class(preds)[2]=="monthly.field.object") {
  i1<-is.element(yy.o,yy.cal)
  i2<-is.element(yy.cal,yy.o)
} else {
  i1<-is.element(yy.o*10000+mm.o*100+dd.o,
                 yy.cal*10000+mm.cal*100+dd.cal)
  i2<-is.element(yy.cal*10000+mm.cal*100+dd.cal,
                 yy.o*10000+mm.o*100+dd.o)
}

# Extract the predictand & predictors:

#print(c(sum(i1),sum(i2)))
y.o<-y.o[i1] ; mm.o<-mm.o[i1]; yy.o<-yy.o[i1]; dd.o<-dd.o[i1]
X.cal<-X.cal[i2,]; mm.cal<-mm.cal[i2]; yy.cal<-yy.cal[i2]; dd.cal<-dd.cal[i2]

# Remove missing values:
i3 <- is.finite(y.o)
#print(c(length(y.o),length(yy.o),length(X.cal[,1]),length(yy.cal)))
y.o<-y.o[i3]; mm.o<-mm.o[i3]; yy.o<-yy.o[i3]; dd.o<-dd.o[i3]
X.cal<-X.cal[i3,]; mm.cal<-mm.cal[i3]; yy.cal<-yy.cal[i3]; dd.cal<-dd.cal[i3]

if (!silent) print("Common times:")
if (!silent) print(range(yy.o))
if (!silent) print(range(y.o))

#--------------------------------------------------------
# De-trend the data used for model calibration:

if (ldetrnd) {
  for (i in 1:length(preds$var.eof)) {
    trnd<-seq(-1,1,length=length(X.cal[,i]))
    dtrnd<-lm(X.cal[,i] ~trnd)
    X.cal[,i]<-dtrnd$residual   
  }
}
if (method=="anm") y <- y.o   # Analog model
              else y <- y.o - mean(y.o,na.rm=TRUE)
trnd<-seq(-1,1,length=length(y))
dtrnd<-lm(y ~ trnd)
if (ldetrnd) {
  y<-dtrnd$residual
}


# Stepwise regression
#scen.gcm.str <- "data.frame("
#calibrate.str <- "data.frame(y=y,"
#for (ipre in 1:length(preds$var.eof)) {
# if (weight) scen.gcm.str <-
#paste(scen.gcm.str,"X",ipre,"=X.gcm[,",ipre,
#                                "]* preds$W[",ipre,"],",sep="")
#else scen.gcm.str <-
#paste(scen.gcm.str,"X",ipre,"=X.gcm[,",ipre,"],",sep="")
 
# Stepwise regression
scen.gcm.str <- "data.frame("
calibrate.str <- "data.frame(y=y,"
for (ipre in 1:length(preds$var.eof)) {
  scen.gcm.str <- paste(scen.gcm.str,"X",ipre,"=X.gcm[,",ipre,
                        "]* preds$W[",ipre,"],",sep="")
  if (method=="anm") {   # Analog model
    if (anm.weight) calibrate.str <- paste(calibrate.str,"X",ipre,"=X.cal[,",
                           ipre,"]* preds$W[",ipre,"],",sep="") else
                    calibrate.str <- paste(calibrate.str,"X",ipre,"=X.cal[,",
                           ipre,"],",sep="") 
  } else calibrate.str <- paste(calibrate.str,"X",ipre,"=X.cal[,",ipre,
                           "]* preds$W[",ipre,"],",sep="")
}
scen.gcm.str <- paste(scen.gcm.str,"yy=yy.gcm,mm=mm.gcm,dd=dd.gcm)",sep="")
scen.gcm <- eval(parse(text=scen.gcm.str))


calibrate.str <- paste(calibrate.str,"yy=yy.cal,mm=mm.cal,dd=dd.cal)",sep="")
calibrate <- eval(parse(text=calibrate.str))

#print(summary(calibrate))
# Due to a bug in step, 'attatch' cannot be used, so it's done
# in a more complicated way.
attach(calibrate)
exprn <- paste(method,"(y ~ 1",sep="")
for (i.eof in 1:length(i.eofs)) {  
  eval(parse(text=
             paste("X",i.eofs[i.eof]," <- calibrate$X",i.eofs[i.eof],sep="")))
}
for (i.eof in 1:length(i.eofs)) {  
  exprn <- paste(exprn," + X",i.eofs[i.eof],sep="")
}
if (method!="anm") exprn <- paste(exprn,xtr.args,")",sep="") else 
                   exprn <- paste(exprn,",","data=calibrate",xtr.args,")",sep="") 

if (!silent) print(paste("Model: ",exprn))
expr <- parse(text=exprn)
lm.mod <- eval(expr)
meths <- methods(class(lm.mod))
if (!silent) print("Stepwise..")
if ((swsm!="none") & !is.null(swsm)) {
  step.wise <- eval(parse(text=paste(swsm,"(lm.mod,trace=0)",sep="")))
}  else step.wise<-lm.mod

#print("ANOVA from step-wise regression:")

stat <- summary(step.wise)
if (length(step.wise$coefficients)>1) {
  if (!is.null(stat$r.squared)) {
    r2 <- round(stat$r.squared*100)
    p.val <- round(100*(1-pf(stat$fstatistic[1],
                           stat$fstatistic[2],
                           stat$fstatistic[3])))
  } else if (method=="anm") {
    cor.test(y,eval(parse(text=paste(predm,"(lm.mod)",sep=""))))
    r2.stat <- cor.test(y,predict.anm(lm.mod))
    r2 <- round(100*r2.stat$estimate^2,2)
    p.val <- round(100*r2.stat$p.value,2)
  } else {
    r2.stat <- eval(parse(text=paste("r2.stat <- cor.test(y,",predm,"(lm.mod))",sep="")))
    r2 <- round(100*r2.stat$estimate^2,2)
    p.val <- round(100*r2.stat$p.value,2)
  }
  fit.p<-as.character(p.val)
} else {
  if (!silent) print("-----------Step failed:----------")
  if (!silent) print(paste("---------",method,":"))
  if (!silent) print(summary(lm.mod))
  if (!silent) print("-----------Step:")
  if (!silent) print(stat)
  r2 <- 0
  p.val <- 100
  fit.p <- "100"
}

# Downscale predictions

#pre.y  <-predict(step.wise)
pre.y  <- eval(parse(text=paste(predm,"(step.wise)",sep="")))
for (i.eof in 1:20) {
  eval(parse(text=paste("rm (X",i.eofs[i.eof],")",sep="")))
}
detach(calibrate)
attach(scen.gcm)
#pre.gcm<-predict(step.wise,newdata=scen.gcm)
pre.gcm <- eval(parse(text=paste(predm,"(step.wise,newdata=scen.gcm)",sep="")))
detach(scen.gcm)
if (!silent) print(summary(step.wise))

# A "fudge" to avoid problems when stepwise rejects all the predictors
# (i.e. only returns an intercept)
if (length(pre.gcm)==1) {
  if (!silent) print(c(length(pre.gcm),c(length(yy.gcm))))
  pre.gcm <- rep(pre.gcm,length(yy.gcm))
}

#print(summary(pre.y))
#print(summary(pre.gcm))
#print(c(mean(pre.gcm[yy.gcm<2010],na.rm=TRUE),
#        mean(y.o[yy.o>1980],na.rm=TRUE)))

if ((class(dat)[2]=="daily.station.record") & (rmac)) {
#  tim.cal <- julian(mm.cal,dd.cal,yy.cal,origin.=c(1,1,1970))
#  tim.cal <- mdy.date(mm.cal, dd.cal, yy.cal)
  tim.cal <- julday(mm.cal, dd.cal, yy.cal) - julday(1,1,1970)
  rm(ac.mod)
  nt.cal <- length(tim.cal)
  ac.mod<-matrix(rep(NA,nt.cal*6),nt.cal,6)
  ac.mod[,1]<-cos(2*pi*tim.cal/365.25); ac.mod[,2]<-sin(2*pi*tim.cal/365.25)
  ac.mod[,3]<-cos(4*pi*tim.cal/365.25); ac.mod[,4]<-sin(4*pi*tim.cal/365.25)
  ac.mod[,5]<-cos(6*pi*tim.cal/365.25); ac.mod[,6]<-sin(6*pi*tim.cal/365.25)
  ac.cal <- data.frame(X=ac.mod)
  rm(ac.mod)
#  tim.gcm <- julian(mm.gcm,dd.gcm,yy.gcm,origin.=c(1,1,1970))
#  tim.gcm <- mdy.date(mm.gcm, dd.gcm, yy.gcm)
  tim.gcm <- julday(mm.gcm, dd.gcm, yy.gcm) - julday(1,1,1970)
  nt.gcm <- length(tim.gcm)
  ac.mod<-matrix(rep(NA,nt.gcm*6),nt.gcm,6)
  ac.mod[,1]<-cos(2*pi*tim.gcm/365.25); ac.mod[,2]<-sin(2*pi*tim.gcm/365.25)
  ac.mod[,3]<-cos(4*pi*tim.gcm/365.25); ac.mod[,4]<-sin(4*pi*tim.gcm/365.25)
  ac.mod[,5]<-cos(6*pi*tim.gcm/365.25); ac.mod[,6]<-sin(6*pi*tim.gcm/365.25)
  ac.gcm <- data.frame(X=ac.mod)
  rm(ac.mod)
#  tim.o <- julian(mm.o,dd.o,yy.o,origin.=c(1,1,1970))
#  tim.o <- mdy.date(mm.o, dd.o, yy.o)
  tim.o <- julday(mm.o, dd.o, yy.o) - julday(1,1,1970)
  nt <- length(tim.o)
  ac.mod<-matrix(rep(NA,nt*6),nt,6)
  ac.mod[,1]<-cos(2*pi*tim.o/365.25); ac.mod[,2]<-sin(2*pi*tim.o/365.25)
  ac.mod[,3]<-cos(4*pi*tim.o/365.25); ac.mod[,4]<-sin(4*pi*tim.o/365.25)
  ac.mod[,5]<-cos(6*pi*tim.o/365.25); ac.mod[,6]<-sin(6*pi*tim.o/365.25)
  ac.obs <- data.frame(X=ac.mod)
  y.o <- y.o + predict(ac.fit,newdata=ac.obs)
  pre.y <- pre.y + predict(ac.fit,newdata=ac.cal)
  pre.gcm <- pre.gcm + predict(ac.fit,newdata=ac.gcm)  
}

if (method!="anm") {
  pre.y   <- pre.y   - mean(pre.y,na.rm=TRUE) + mean(y.o,na.rm=TRUE)
  pre.gcm <- pre.gcm - mean(pre.gcm[yy.gcm<2010],na.rm=TRUE) +
                     mean(y.o[yy.o>1980],na.rm=TRUE)
}

if (!silent) print(summary(pre.gcm))

if ( (regexpr("precip",lower.case(dat$obs.name)) > 0) |
     (regexpr("rain",lower.case(dat$obs.name)) > 0) ) {
  pre.y[pre.y < 0] <-  0
  pre.gcm[pre.gcm < 0] <-  0
}

# Predictions: GCM
# Determine which PCs were selected in the step procedure
# Note, some of the higher modes are truncated

c<-as.character(step.wise$call[2])
c<-unlist(strsplit(c," \\~ "))
c<-c[2]
c<-paste(unlist(strsplit(c," \\+ ")),' ')
#print(c)
incl<- rep(FALSE,length(preds$var.eof))
for (i in 1:length(i.eofs)) {
  if (!is.na( charmatch(paste('X',as.character(i),' ',sep=""),c ) ))  {
    incl[i]<-TRUE
  }
}
#print(incl)

# Note that the intercept is included in lm.coe, but not in the
# coefficients held by c.

lm.coe <- coef(step.wise)

# Find the predictor patterns

preds2D<-preds$EOF
dims <- dim(preds2D)
if (length(dims) > 2) dim(preds2D)<-c(dims[1],dims[2]*dims[3])

if (!silent) print("Reconstruct the spatial patterns")
if (!silent) print(lm.coe)

i.last <- 0
list.expr <- "list("
id <- row.names(table(preds$id.x))
#print(id)
#print(table(preds$id.lon))
#print(table(preds$id.lat))
for (i in 1:n.fld) {
#  print(id[i])
  i.lon <- preds$id.lon == id[i]
  i.lat <- preds$id.lat == id[i]
  ny<-preds$size[2,i]
  nx<-preds$size[3,i]
  i.fld <- seq(i.last+1,i.last+ny*nx,by=1)
  i.last <- max(i.fld)
#  print(paste("Dimension of field ",i))
#  print(dim(preds2D))
#  print(c(sum(incl),sum(i.fld)))
  EOF.1 <- t(preds2D[,i.fld])
  EOF.1 <-  EOF.1[,incl]
#  print(dim(EOF.1))
  expr <- paste("X.",i," <- cbind(0,EOF.1) %*% lm.coe[1:(sum(incl)+1)]",sep="")
  eval(parse(text=expr))
#  print(paste("2D -> 3D: nx=",nx," ny=",ny))
#  print(eval(parse(text=paste("dim(X.",i,")",sep=""))))
  expr <- paste("dim(X.",i,") <- c(ny,nx)",sep="")
  eval(parse(text=expr))
  eval(parse(text=paste("lon.",i," <- preds$lon[i.lon]",sep="")))
  eval(parse(text=paste("lat.",i," <- preds$lat[i.lat]",sep="")))
  list.expr <- paste(list.expr,"X.",i,"=X.",i,
                     ", lon.",i,"=lon.",i,
                     ", lat.",i,"=lat.",i,", ",sep="")
#  print(length(lon.1))
#  print(length(lat.1))
}

# Linear trend:

if (!silent) print("Linear trend for GCM (deg C/decade)")
x.ind <- seq(0,1,length=length(yy.gcm))
tr.dat<-data.frame(y=pre.gcm, x=x.ind)
nt <- length(yy.gcm)
lm.tr <- lm(y ~ x, data=tr.dat)
stat.tr.fit <- summary(lm.tr)
coef.fit<-stat.tr.fit$coefficients
rate.ds <- round(as.real(round(coef.fit[2]*10,2))*(x.ind[2]-x.ind[1]),2)
rate.err  <- round(as.real(round(coef.fit[4]*10,2))*(x.ind[2]-x.ind[1]),2)

#print(coef.fit)
if (!silent) print("Slope and its uncertainty")
if (!silent) print(c(rate.ds,rate.err))
pre.fit<-predict(lm.tr,data= yy)

# Polinomial trend

#print("Polinomial trend")
lm.tr.p<-lm(y ~ x + I(x^2) +I(x^3) + I(x^4) + I(x^5), data=tr.dat)
pre.p.fit<-predict(lm.tr.p,data=tr.dat)
coef.p.fit<-lm.tr.p$coefficients
coef.p.fit[is.na(coef.p.fit)] <- 0
der.p.fit<-c(coef.p.fit[2],2*coef.p.fit[3],3*coef.p.fit[4],
             4*coef.p.fit[5],5*coef.p.fit[6])*(x.ind[2]-x.ind[1])
tr.est.p.fit<-(der.p.fit[1] + der.p.fit[2]*x.ind + der.p.fit[3]*x.ind^2 +
               der.p.fit[4]*x.ind^3 + der.p.fit[5]*x.ind^4)*10
gcm.stat <- summary(lm.tr)

# Estimate the P-values associated with the trends

if (!is.null(gcm.stat$fstatistic)) {
  gcm.trnd.p<-as.character(round(100*(1-pf(gcm.stat$fstatistic[1],
                                           gcm.stat$fstatistic[2],
                                           gcm.stat$fstatistic[3]))))
} else {gcm.trnd.p<-"100"}
gcm.trnd.r2 <- gcm.stat$r.squared

if (!silent) print(paste("P-value of fit=",fit.p))
if (!silent) print(paste("P-value of trend-fit for downscaled scenario",gcm.trnd.p))

#---------------------------------------------------

if ((method!="nnet") & (method!= "anm") & lsave) {
  
mod.name<-paste(direc,"ds.mod_",preds.id,"_",region,"_",
             substr(dat$location,1,eos),"_",dat$ele,"_",preds$c.mon,'_',
             substr(attr(preds$tim,"unit"),1,3),"_",method,
             ex.tag,sep="")

mod.tab <- xtable(step.wise,
                  caption=paste("Calibration period: ",month,
                  " ",range(yy.cal)[1],"-",range(yy.cal)[2],
                    " using ",preds$id.t[cal.id],sep=""))
print.xtable(mod.tab,type="latex",
             file=paste(mod.name,".tex",sep=""))
print.xtable(mod.tab,type="html",
             file=paste(mod.name,".html",sep=""))
}
sce.name<-paste(direc,"ds.res_",preds.id,"_",region,"_",
             substr(dat$location,1,eos),"_",dat$ele,"_",preds$c.mon,'_',
             substr(attr(preds$tim,"unit"),1,3),"_",method,
             ex.tag,sep="")

scen.table<-xtable(data.frame(year=yy.gcm,
                              downscaled=round(pre.gcm,2)),
                   caption=paste("Linear trend=",rate.ds,
                     ds.unit,"/decade over ",month," ",
                     range(yy.gcm)[1],"-",range(yy.gcm)[2],
                     " using",preds$id.t[preds$id.t!=cal.id][1],
                     "; p-value for linear trend-fit=",
                     gcm.trnd.p,"%.",sep=""))

if (lsave) print.xtable(scen.table,type="html",
           file=paste(sce.name,".html",sep=""))

pred.name <- row.names(table(preds$id.x))
list.expr <- paste(list.expr,
         "lon.loc=dat$lon,lat.loc=dat$lat,alt.loc=dat$alt,",
         "step.wise=step.wise,location=loc,",
         "yy.gcm=yy.gcm, mm.gcm=mm.gcm, dd.gcm=dd.gcm, ",
         "yy.cal=yy.cal, dd.cal=dd.cal,","mm.cal=mm.cal,",
         "n.fld=n.fld,unit=ds.unit,",
         "rate.ds=rate.ds,rate.err=rate.err,gcm.trnd.p=gcm.trnd.p,",
         "y.o=y.o,mm.o=mm.o,yy.o=yy.o,dd.o=dd.o,",
         "fit.p=fit.p,fit.r2=r2,pre.p.fit=pre.p.fit,",
         "pre.gcm=pre.gcm,pre.y=pre.y,gcm.stat=gcm.stat,",
         "month=month,v.name=v.name, region=preds$region,",
         "id.1=cal.id,id.2=preds$id.t[preds$id.t!=cal.id][1],",
         "pre.fit=pre.fit,tr.est.p.fit=tr.est.p.fit,ex.tag=ex.tag,",
         "pred.name=pred.name)",sep="")       
#print(list.expr)
ds<-eval(parse(text=list.expr))
if (!silent) print(paste("File name:",fname))
class(ds) <- "ds"
if (lsave) save(file=fname,ds,ascii=FALSE) 
#print("Plotting...")
#print(preds$region)
if (plot) plotDS(ds,leps)
invisible(ds)
}
