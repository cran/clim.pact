rm(list=ls())
library(clim.pact)
source("clim.pact/R/eof.R")
source("clim.pact/R/ds.R")
source("clim.pact/R/cat.fields.R")
source("clim.pact/R/mix.fields.R")

cmon<-c("Jan","Feb","Mar","Apr","May","Jun",
        "Jul","Aug","Sep","Oct","Nov","Dec")

domains <- c("D1","D2","D3","D4","D5","D6","D7","D8")
analyses  <- c("DNMI","NCEP")

for (domain in domains) {
  for (analysis in analyses) {

x.domain <- switch(domain,
                "D1"=list(x.rng=c(-90,60),y.rng=c(30,80)),
                "D2"=list(x.rng=c(-60,40),y.rng=c(40,70)),
                "D3"=list(x.rng=c(-40,40),y.rng=c(50,70)),
                "D4"=list(x.rng=c(-20,40),y.rng=c(50,70)),
                "D5"=list(x.rng=c(  0,35),y.rng=c(55,70)),
                "D6"=list(x.rng=c( 10,50),y.rng=c(60,75)),
                "D7"=list(x.rng=c(-20,20),y.rng=c(55,65)),
                "D8"=list(x.rng=c(-40,20),y.rng=c(65,85)))
x.rng <- x.domain$x.rng
y.rng <- x.domain$y.rng

############################## DNMI ###########################

if (analysis=="DNMI") {
  print("###------------ DNMI -------------###")
  t.rng <- c(1900,2099)
  obs.t2m <- retrieve.nc("/home/kareb/data/analysis/T2M_p.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=t.rng)
  obs.slp <- retrieve.nc("/home/kareb/data/analysis/SLP_oi.nc",
                         x.rng=x.rng,y.rng=y.rng,t.rng=t.rng)
} else {

############################## NCEP ###########################

  print("###------------ NCEP -------------###")
  t.rng <- c(1958,2099)
  obs.t2m <- retrieve.nc("/home/kareb/data/ncep/ncep_t2m.nc",
                         x.rng=x.rng,y.rng=y.rng,t.rng=t.rng)
  obs.slp <- retrieve.nc("/home/kareb/data/ncep/ncep_slp.nc",
                         x.rng=x.rng,y.rng=y.rng,t.rng=t.rng)
}

##################### HadCM3: B2 ########################

print("### HadCM3 B2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/HADCM3_B2_tem.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=t.rng)
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/HADCM3_B2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=t.rng)

for (im in 1:12) {
  print("cat.fields")
  XX.t2m <- cat.fields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- eof(XX.t2m,plot=FALSE)
  XX.slp <- cat.fields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- eof(XX.slp,plot=FALSE)
  print("mix.fields:")
  XX <- mix.fields(XX.t2m,XX.slp)
  print("mix.fields EOF")
  EOF.t2m.slp <- eof(XX,plot=FALSE)
}

##################### HadCM3: A2 ########################

print("### HadCM3 A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/HADCM3_A2_tem.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2000,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/HADCM3_A2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2000,2099))

for (im in 1:12) {
  print("cat.fields")
  XX.t2m <- cat.fields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- eof(XX.t2m,plot=FALSE)
  XX.slp <- cat.fields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- eof(XX.slp,plot=FALSE)
  print("mix.fields:")
  XX <- mix.fields(XX.t2m,XX.slp)
  print("mix.fields EOF")
  EOF.t2m.slp <- eof(XX,plot=FALSE)
}

##################### ECHAM4/OPYC3: B2 ########################

print("### ECHAM4/OPYC3 B2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/EH4OPYC_B2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/EH4OPYC_B2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
for (im in 1:12) {
  XX.t2m <- cat.fields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- eof(XX.t2m,plot=FALSE)
  XX.slp <- cat.fields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- eof(XX.slp,plot=FALSE)
  XX <- mix.fields(XX.t2m,XX.slp)
  EOF.t2m.slp <- eof(XX,plot=FALSE)
}

##################### ECHAM4/OPYC3: A2 ########################

print("### ECHAM4/OPYC3 A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/EH4OPYC_A2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/EH4OPYC_A2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
for (im in 1:12) {
  XX.t2m <- cat.fields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- eof(XX.t2m,plot=FALSE)
  XX.slp <- cat.fields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- eof(XX.slp,plot=FALSE)
  XX <- mix.fields(XX.t2m,XX.slp)
  EOF.t2m.slp <- eof(XX,plot=FALSE)
}

##################### CSIRO: B2 ########################


print("### CSIRO B2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/CSIRO_B2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1961,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/CSIRO_B2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1961,2099))
for (im in 1:12) {
  XX.t2m <- cat.fields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- eof(XX.t2m,plot=FALSE)
  XX.slp <- cat.fields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- eof(XX.slp,plot=FALSE)
  XX <- mix.fields(XX.t2m,XX.slp)
  EOF.t2m.slp <- eof(XX,plot=FALSE)
}

##################### CSIRO: A2 ########################

print("### CSIRO A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/CSIRO_A2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1961,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/CSIRO_A2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1961,2099))
for (im in 1:12) {
  XX.t2m <- cat.fields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- eof(XX.t2m,plot=FALSE)
  XX.slp <- cat.fields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- eof(XX.slp,plot=FALSE)
  XX <- mix.fields(XX.t2m,XX.slp)
  EOF.t2m.slp <- eof(XX,plot=FALSE)
}


##################### CCC: B2 ########################

print("### CCC B2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/CCCma_B2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/CCCma_B2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
for (im in 1:12) {
  XX.t2m <- cat.fields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- eof(XX.t2m,plot=FALSE)
  XX.slp <- cat.fields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- eof(XX.slp,plot=FALSE)
  XX <- mix.fields(XX.t2m,XX.slp)
  EOF.t2m.slp <- eof(XX,plot=FALSE)
}

##################### CCC: A2 ########################

print("### CCC A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/CCCma_A2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/CCCma_A2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
for (im in 1:12) {
  XX.t2m <- cat.fields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- eof(XX.t2m,plot=FALSE)
  XX.slp <- cat.fields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- eof(XX.slp,plot=FALSE)
  XX <- mix.fields(XX.t2m,XX.slp)
  EOF.t2m.slp <- eof(XX,plot=FALSE)
}


##################### NCARCSM: A2 ########################

print("### NCARCSM: A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/NCARCSM_A2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/NCARCSM_A2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
for (im in 1:12) {
  XX.t2m <- cat.fields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- eof(XX.t2m,plot=FALSE)
  XX.slp <- cat.fields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- eof(XX.slp,plot=FALSE)
  XX <- mix.fields(XX.t2m,XX.slp)
  EOF.t2m.slp <- eof(XX,plot=FALSE)
}


##################### NCARPCM: B2 ########################

print("### NCARPCM: B2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/NCARPCM_B2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/NCARPCM_B2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
for (im in 1:12) {
  XX.t2m <- cat.fields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- eof(XX.t2m,plot=FALSE)
  XX.slp <- cat.fields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- eof(XX.slp,plot=FALSE)
  XX <- mix.fields(XX.t2m,XX.slp)
  EOF.t2m.slp <- eof(XX,plot=FALSE)
}

##################### NCARPCM: A2 ########################


print("### NCARPCM: A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/NCARPCM_A2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/NCARPCM_A2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
for (im in 1:12) {
  XX.t2m <- cat.fields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- eof(XX.t2m,plot=FALSE)
  XX.slp <- cat.fields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- eof(XX.slp,plot=FALSE)
  XX <- mix.fields(XX.t2m,XX.slp)
  EOF.t2m.slp <- eof(XX,plot=FALSE)
}

##################### GFDL: B2 ########################

print("### GFDL: A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/GFDL_B2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/GFDL_B2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
for (im in 1:12) {
  XX.t2m <- cat.fields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- eof(XX.t2m,plot=FALSE)
#  XX.slp <- cat.fields(obs.slp,gcm.slp,mon=im)
#  EOF.slp <- eof(XX.slp,plot=FALSE)
#  XX <- mix.fields(XX.t2m,XX.slp)
#  EOF.t2m.slp <- eof(XX,plot=FALSE)
}

##################### GFDL: A2 ########################

print("### GFDL: A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/GFDL_A2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/GFDL_A2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
for (im in 1:12) {
  XX.t2m <- cat.fields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- eof(XX.t2m,plot=FALSE)
#  XX.slp <- cat.fields(obs.slp,gcm.slp,mon=im)
#  EOF.slp <- eof(XX.slp,plot=FALSE)
#  XX <- mix.fields(XX.t2m,XX.slp)
#  EOF.t2m.slp <- eof(XX,plot=FALSE)
}

}
}

