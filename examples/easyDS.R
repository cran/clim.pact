rm(list=ls())
library(clim.pact)
#source("clim.pact/R/getnordklim.R")

cmon<-c("Jan","Feb","Mar","Apr","May","Jun",
        "Jul","Aug","Sep","Oct","Nov","Dec")
ele.c <- c("TAM","TAX","Th","Thd","TAN","Tl","Tld","SLP","RR","RRX",
           "DSC","CLOUD","SDM")
ele <- c(101,111,112,113,121,122,123,401,601,602,701,801,911)

if (!file.exists("data/avail.nordklim.Rdata")) {
  print("Please be patient - searching through the files to check data")
  locs <- avail.locs()
  names <- locs$name[locs$ident=="NORDKLIM"]
  names <- names[!is.na(names)]
  x <- rep(FALSE,length(names))
  evalstr <- "list(name=names"
  for (i in 1:length(ele)) evalstr <- paste(evalstr,", ",ele.c[i],"= x",sep="")
  evalstr <- paste(evalstr,")",sep="")
  print(evalstr)
  avail.nordklim <- eval(parse(text=evalstr))
  for (iele in 1:length(ele)) {
    for (iloc in 1:length(names)) {
      obs <- getnordklim(names[iloc],ele.c=ele[iele])
      if (is.finite(min(obs$yy))) {
        eval(parse(text=paste("avail.nordklim$",ele.c[iele],"[iloc] <- T",sep="")))
        print(c(names[iloc],ele.c[iele],min(obs$yy),max(obs$yy)))
      }
    }
  }
  save(file="data/avail.nordklim.Rdata",avail.nordklim)
} else load("data/avail.nordklim.Rdata")
                                     
sce <- c("A1","B1","A2","B2")
print(sce)
sce <- sce[as.integer(readline(prompt="Which SRES scenario (#)?"))]
elems <- avail.elem()
print(elems$name)
iele <- as.integer(readline(prompt="Which element (#)?"))

ii <- eval(parse(text=paste("avail.nordklim$",ele.c[iele],sep="")))
locs <- avail.nordklim$name[ii]
ident <- rep("NORDKLIM",length(locs))
if ((iele==1) | (iele==9)) {
  all.locs <- avail.locs()
  locs <- c(locs,all.locs$name[all.locs$ident=="NACD"])
  ident <- c(ident,all.locs$ident[all.locs$ident=="NACD"])
}
print(locs)
iloc <- as.integer(readline(prompt="Which location (#)?"))

mon <- as.integer(readline(prompt="Which month (1-12)?"))
preds <- avail.preds()
preds <- preds[grep(cmon[mon],preds)]
preds <- preds[grep(sce,preds)]
print(preds)
ipre <- as.integer(readline(prompt="Which predictor (#)?"))

print("You selected:")
print(paste(locs[iloc]," & ",preds[ipre]," & ",ele.c[iele]," & ",ident[iloc]))
load(paste("data/",preds[ipre],sep=""))
if (ident[iloc]=="NACD") obs <- getnacd(locs[iloc],ele.c=as.character(ele[iele])) else
                         obs <- getnordklim(locs[iloc],ele.c=as.character(ele[iele]))
ds <- DS(preds=eof,dat=obs,plot=FALSE)
plot(range(c(ds$yy.gcm,ds$yy.cal)),type="n",
     range(c(ds$pre.gcm,ds$pre.y)),
     sub=paste("Calibration: ",ds$id.1,", Scenario: ",ds$id.2,sep=""),
     main=obs$location,xlab="Time",
     ylab=paste(obs$obs.name,"   [",obs$unit,"]",sep=""))
grid()
ii <- ds$yy.gcm > max(ds$yy.cal)
lines(ds$yy.gcm[ii],ds$pre.gcm[ii],type="s",lwd=3,col="darkred")
lines(ds$yy.gcm,ds$pre.gcm,type="s",lwd=1,lty=2,col="darkred")
lines(ds$yy.cal,ds$pre.y,type="s",lwd=3,col="darkblue")
lines(ds$yy.gcm,ds$pre.fit, col = "red",lwd=1,lty=2)

text(quantile(c(ds$yy.cal,ds$yy.gcm),0.01),
     min(c(ds$pre.y,ds$pre.gcm)),pos=4,cex=0.6,
     paste(ds$month,": Trend fit: P-value=",ds$gcm.trnd.p,"%; ",
           "Projected trend= ",ds$rate.ds,"+-",ds$rate.err," ",
           ds$unit,"/decade",sep=""))


