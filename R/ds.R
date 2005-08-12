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
               swsm="step",predm="predict",lsave=FALSE,rmac=TRUE,
               silent=FALSE) {


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

#REB 09.03.05
#print(summary(c(dat$val)))
#print(paste(sum(is.na(dat$val)),"points are allready NA"))
#print(paste("Setting",sum(dat$val < -99,na.rm=TRUE),"points to NA"))

dat$val[dat$val < -99] <-NA

if (rmac) {
   #print("TEST: remove the annual cycle with anomaly.station")
   dat <- anomaly.station(dat)
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
  if (!silent) print(attr(preds$tim,"unit"))
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
#  v.name  <- abbreviate(dat$obs.name)
  v.name  <- dat$obs.name
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

  if (!is.null(preds$attributes$daysayear)) daysayear <- preds$attributes$daysayear else
                                            daysayear <- 365.25

if (class(dat)[2]=="monthly.station.record"){
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
}

#mm.o <- mm.o[!is.na(y.o)]
#yy.o <- yy.o[!is.na(y.o)]
#dd.o <- dd.o[!is.na(y.o)]
#y.o <- y.o[!is.na(y.o)]
#REB 09.03.05
#print(paste(">---1: length y.o=",length(y.o),"length(yy.o)=",length(yy.o)))

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
#print(paste(">---2: length y.o=",length(y.o),"length(yy.o)=",length(yy.o)))

if (is.null(cal.id)) cal.id<- preds$id.t[1]

#REB 09.03.05
#print(dim(preds$PC))
#print(length(preds$id.t))
#print(paste("Calibration predictors:",cal.id))
#print(sum(preds$id.t==cal.id & !is.na(preds$PC[,1])))
#print(range(preds$yy[preds$id.t==cal.id & !is.na(preds$PC[,1])]))

X.cal<-  preds$PC[preds$id.t==cal.id & !is.na(preds$PC[,1]),]
yy.cal<- preds$yy[preds$id.t==cal.id & !is.na(preds$PC[,1])]
mm.cal<- preds$mm[preds$id.t==cal.id & !is.na(preds$PC[,1])]
dd.cal<- preds$dd[preds$id.t==cal.id & !is.na(preds$PC[,1])]

if (!silent) print("------------Match times---------------- ")
#REB 09.03.05
#print(range(yy.cal))
X.cal<-  X.cal[is.element(mm.cal,mon),]
yy.cal<- yy.cal[is.element(mm.cal,mon)]
dd.cal<- dd.cal[is.element(mm.cal,mon)]
mm.cal<- mm.cal[is.element(mm.cal,mon)]
#print(c(length(y.o),length(yy.o),length(X.cal[,1]),length(yy.cal)))

if (sum((preds$id.t!=cal.id) & !is.na(preds$PC[,1]))>0) {
  X.gcm<-preds$PC[preds$id.t!=cal.id  & !is.na(preds$PC[,1]),]
  yy.gcm<-as.vector(preds$yy[preds$id.t!=cal.id & !is.na(preds$PC[,1])])
  mm.gcm<-as.vector(preds$mm[preds$id.t!=cal.id & !is.na(preds$PC[,1])])
  dd.gcm<-as.vector(preds$dd[preds$id.t!=cal.id & !is.na(preds$PC[,1])])
  
#REB 09.03.05
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
#REB 09.03.05
#print(paste(">---3: length y.o=",length(y.o),"length(yy.o)=",length(yy.o)))

if ((class(preds)[2]=="monthly.field.object") |
    (class(preds)[3]=="monthly.field.object")) {
  i1<-is.element(yy.o,yy.cal)
  i2<-is.element(yy.cal,yy.o)
} else {
  i1<-is.element(yy.o*10000+mm.o*100+dd.o,
                 yy.cal*10000+mm.cal*100+dd.cal)
  i2<-is.element(yy.cal*10000+mm.cal*100+dd.cal,
                 yy.o*10000+mm.o*100+dd.o)
}

# Extract the predictand & predictors:

if (!silent) print(paste("Number of coinciding obs:",sum(i1),",",sum(i2)))
if (!silent) print(summary(yy.o))
if (!silent) print(summary(yy.cal))

y.o<-y.o[i1] ; mm.o<-mm.o[i1]; yy.o<-yy.o[i1]; dd.o<-dd.o[i1]; y <- y.o
X.cal<-X.cal[i2,]; mm.cal<-mm.cal[i2]; yy.cal<-yy.cal[i2]; dd.cal<-dd.cal[i2]
#print(paste(">---4: length y.o=",length(y.o),"length(y)=",length(y)))


# Remove missing values:
valid.cal <- is.finite(y.o)
#REB 09.03.05
#print(c(length(y.o),length(yy.o),length(X.cal[,1]),length(yy.cal)))
#y<-y.o[valid.cal]; 
#mm.o<-mm.o[valid.cal]; yy.o<-yy.o[valid.cal]; dd.o<-dd.o[valid.cal]
#X.cal<-X.cal[valid.cal,]; mm.cal<-mm.cal[valid.cal]; yy.cal<-yy.cal[valid.cal]; dd.cal<-dd.cal[valid.cal]

if (!silent) print("Common times:")
if (!silent) print(range(yy.o))
if (!silent) {print(summary(y.o))
              print(paste(sum(is.finite(y.o)),"valid points,",sum(!is.finite(y.o)),"NA's"))}

#--------------------------------------------------------
# De-trend the data used for model calibration:

if (ldetrnd) {
  if (!silent) print("de-trend:")
  for (i in 1:length(preds$var.eof)) {
    trnd<-seq(-1,1,length=length(X.cal[,i]))
    dtrnd<-lm(X.cal[,i] ~trnd)
    X.cal[,i]<-dtrnd$residual   
  }
  trnd<-seq(-1,1,length=length(y))
  dtrnd<-lm(y ~ trnd)
  y[!valid.cal] <- NA
  y[valid.cal]<-dtrnd$residual
}

# Assign calibration and prediction data:
n.eofs<- min(c(length(preds$var.eof),length(i.eofs)))    
scen.gcm.str <- "data.frame("
calibrate.str <- "data.frame(y=y,"
valid.cal2 <- rep(TRUE,length(y))
for (ipre in 1:n.eofs) {
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
scen.gcm.str <- paste(scen.gcm.str,"yy=as.vector(yy.gcm),mm=as.vector(mm.gcm),dd=as.vector(dd.gcm))",sep="")
#print("GCM:")
#print(scen.gcm.str)
scen.gcm <- eval(parse(text=scen.gcm.str))

calibrate.str <- paste(calibrate.str,"yy=as.vector(yy.cal),mm=as.vector(mm.cal),dd=as.vector(dd.cal))",sep="")
#print("Calibration:")
#print(calibrate.str); print(c(length(y),NA,dim(X.cal)))
calibrate <- eval(parse(text=calibrate.str))

#print(summary(calibrate))
# Due to a bug in step, 'attatch' cannot be used, so it's done
# in a more complicated way.
attach(calibrate)

# Stepwise regression
if ((!silent) & (method!="anm")) print("stepwise regression:")
exprn <- paste(method,"(y ~ 1",sep="")
for (i.eof in 1:n.eofs) {  
  eval(parse(text=
             paste("X",i.eofs[i.eof]," <- calibrate$X",i.eofs[i.eof],sep="")))
  valid.cal2 <- valid.cal2 & is.finite(eval(parse(text=paste("calibrate$X",i.eofs[i.eof],sep=""))))
}
for (i.eof in 1:n.eofs) {  
  exprn <- paste(exprn," + X",i.eofs[i.eof],sep="")
}
if (method!="anm") exprn <- paste(exprn,xtr.args,",data=calibrate)",sep="") else 
                   exprn <- paste(exprn,",","data=calibrate",xtr.args,")",sep="") 

if (!silent) {print(paste("Model: ",exprn))
              print(paste(sum(valid.cal),"/",sum(valid.cal2),"valid calibr. points"))}

#print(c(length(y),length(X1),length(X2),length(X3),NA,sum(is.finite(calibrate$y)),sum(is.finite(calibrate$X1))))
lm.mod <- eval(parse(text=exprn))
#print(summary(lm.mod))
lm.mod$coefficients[!is.finite(lm.mod$coefficients)] <- 0
meths <- methods(class(lm.mod))
if (!silent) print(paste("Stepwise:   ",swsm,"(lm.mod,trace=0)",sep=""))
if ((swsm!="none") & !is.null(swsm)) {
  step.wise <- eval(parse(text=paste(swsm,"(lm.mod,trace=0)",sep="")))
}  else step.wise<-lm.mod
step.wise$coefficients[!is.finite(step.wise$coefficients)] <- 0

#print(paste(">---5: length y.o=",length(y.o),"length(yy.o)=",length(yy.o)))
if (!silent) print("ANOVA from step-wise regression:")

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
    r2.stat <- eval(parse(text=paste("r2.stat <- cor.test(y,",
                            predm,"(lm.mod))",sep="")))
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
pre.y <- rep(NA,length(yy.cal))
pre.y[valid.cal]  <- eval(parse(text=paste(predm,"(step.wise)",sep="")))
for (i.eof in 1:20) {
  eval(parse(text=paste("rm (X",i.eofs[i.eof],")",sep="")))
}
#print(paste(">---6: length yy.cal=",length(yy.cal),"length(pre.y)=",length(pre.y)))
detach(calibrate)
attach(scen.gcm)
#pre.gcm<-predict(step.wise,newdata=scen.gcm)
if (method=="anm") {           # REB 18.03.2005
  analog <-  eval(parse(text=paste(predm,"(step.wise,se.fit=TRUE)",sep="")))
  pre.gcm <- eval(parse(text=paste(predm,"(step.wise,newdata=scen.gcm)",sep="")))
  analog$yy <- yy.o[analog$date.min]
  analog$mm <- mm.o[analog$date.min]
  analog$dd <- dd.o[analog$date.min]
} else {
  pre.gcm <- eval(parse(text=paste(predm,"(step.wise,newdata=scen.gcm)",sep="")))
}

if (!silent) print("Downscaled anomalies:")
if (!silent) print(summary(pre.gcm))

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

if (rmac) {
  if (!silent) print("Adding the annual cycle to the prediction: preClim.station")
  #print("y.o:")
  y.o <- y.o + preClim.station(dat,dd.o,mm.o,yy.o)
  #print("pre.y:")
  pre.y <- pre.y + preClim.station(dat,dd.cal,mm.cal,yy.cal)
  #print("pre.gcm:")
  pre.gcm <- pre.gcm + preClim.station(dat,dd.gcm,mm.gcm,yy.gcm)
}

if (method!="anm") {
  ii1 <- is.element(yy.gcm,yy.o)
  ii2 <- is.element(yy.o,yy.gcm)
  cal.mean <- mean(pre.y,na.rm=TRUE)
  gcm.mean <- mean(pre.gcm[ii1],na.rm=TRUE)
  obs.mean <- mean(y.o[ii2],na.rm=TRUE)
  obs.mean2 <- mean(y.o,na.rm=TRUE)
  if (!is.finite(gcm.mean)) gcm.mean  <- gcm.mean[1]
  if (!is.finite(obs.mean)) obs.mean  <- 0 
  if (!is.finite(cal.mean)) cal.mean  <- 0
  if (!is.finite(obs.mean2)) obs.mean2  <- obs.mean
  pre.y   <- pre.y   - cal.mean + obs.mean2
#  pre.gcm <- pre.gcm - gcm.mean + obs.mean2   # REB 26.08.2004: 'obs.mean' replaced with 'obs.mean2'
  pre.gcm <- pre.gcm - cal.mean + obs.mean2   # REB 02.02.2005: 'gcm.mean' replaced with 'cal.mean'
  #print(paste(">---7: length y.o=",length(y.o),"length(yy.o)=",length(yy.o)))                     
}
#print("Check#3:")
#print(summary(pre.gcm))  # TEST

#print(paste("DS:dat$obs.name =",dat$obs.name))
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

#print(paste(">---8: length y.o=",length(y.o),"length(yy.o)=",length(yy.o)))
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

#print(paste(">---9: length y.o=",length(y.o),"length(yy.o)=",length(yy.o)))

if (!silent) print("Dignosis:")
if (!silent) print(direc)
if (!silent) print(preds.id)
if (!silent) print(region)
if (!silent) print(dat$location)
if (!silent) print(dat$ele)
if (!silent) print(preds$c.mon)
if (!silent) print(attr(preds$tim,"unit"))
if (!silent) print(method)


sce.name<-paste(direc,"ds.res_",preds.id,"_",region,"_",
             substr(dat$location,1,eos),"_",dat$ele,"_",preds$c.mon,'_',
             substr(attr(preds$tim,"unit"),1,3),"_",method,
             ex.tag,sep="")


#print(paste(">---10: length y.o=",length(y.o),"length(yy.o)=",length(yy.o)))

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
         "pred.name=pred.name,sce.name=sce.name,preds.name=preds$f.name)",sep="")       
#print(list.expr)
ds<-eval(parse(text=list.expr))
if (!silent) print(paste("File name:",fname))
if (method=="anm") ds$analog <- analog
class(ds) <- "ds"
if (lsave) save(file=fname,ds,ascii=FALSE) 
#print("Plotting...")
#print(preds$region)
if (plot) {
   ds.map <- plotDS(ds,leps)
   ds$map <- ds.map
} 
invisible(ds)
}
