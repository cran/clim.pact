# NCAR PCM A2 :      1980-2099 (Jan-Feb 1980 empty)
# NCAR PCM B2 :      2000-2099
# NCAR CSM A2 :      2000-2099
# CSIRO A1/B1/A2/B2: 1961-2100
# CCCma A2/B2:       1900-2100
# ECHAM4 A2/B2:      1990-2100
# HadCM3 A2/B2:      1950-2099


rm(list=ls())
library(clim.pact)

do.all <- FALSE
cmon<-c("Jan","Feb","Mar","Apr","May","Jun",
        "Jul","Aug","Sep","Oct","Nov","Dec")

#domains <- c("D2","D4","D9","D10","D1","D3","D5","D6","D7","D8","D11")
#analyses  <- c("NCEP","DNMI")
domains <- c("D14","D4","D17","D5","D11","D16")
analyses  <- c("NCEP")

for (domain in domains) {
  for (analysis in analyses) {

x.domain <- switch(domain,
                "D1"=list(x.rng=c(-90,60),y.rng=c(30,80)),
                "D2"=list(x.rng=c(-60,40),y.rng=c(40,70)),
                "D3"=list(x.rng=c(-40,40),y.rng=c(50,70)),
                "D4"=list(x.rng=c(-20,40),y.rng=c(50,70)),
                "D5"=list(x.rng=c(  0,35),y.rng=c(55,70)),
                "D6"=list(x.rng=c( -10,50),y.rng=c(50,75)),
                "D7"=list(x.rng=c(-20,20),y.rng=c(55,65)),
                "D8"=list(x.rng=c(-40,20),y.rng=c(65,85)),
                "D9"=list(x.rng=c(-60,-40),y.rng=c(40,70)),
                "D10"=list(x.rng=c(-90,-30),y.rng=c(60,80)),
                "D11"=list(x.rng=c(-40,40),y.rng=c(50,80)),
                "D12"=list(x.rng=c(-70,0),y.rng=c(60,80)),
                "D13"=list(x.rng=c(-90,-30),y.rng=c(50,80)),
                "D14"=list(x.rng=c(0,50),y.rng=c(60,83)),
                "D15"=list(x.rng=c(-20,70),y.rng=c(50,83)),
                "D16"=list(x.rng=c(0,50),y.rng=c(50,83)),
                "D17"=list(x.rng=c(-20,50),y.rng=c(30,75)),
                "D18"=list(x.rng=c(-35,40),y.rng=c(65,85)))
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

##################### HadCM3: B2 ##############################

print("### HadCM3 B2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/HADCM3_B2_tem.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1950,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/HADCM3_B2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1950,2099))

if (do.all) {  # REB 30.09.03

for (im in 1:12) {
  print("catFields")
  XX.t2m <- catFields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- EOF(XX.t2m,plot=FALSE)
  XX.slp <- catFields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- EOF(XX.slp,plot=FALSE)
  print("mixFields:")
  XX <- mixFields(XX.t2m,XX.slp)
  print("mixFields EOF")
  EOF.t2m.slp <- EOF(XX,plot=FALSE)
}

} # end do REB 30.09.03

gcm.t2m.0 <- gcm.t2m
gcm.slp.0 <- gcm.slp

##################### HadCM3: A2 ########################

print("### HadCM3 A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/HADCM3_A2_tem.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2000,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/HADCM3_A2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2000,2099))

gcm.t2m <- catFields(gcm.t2m.0,gcm.t2m,
                     interval.1=c(1950,1999),interval.2=c(2000,2099))
gcm.slp <- catFields(gcm.slp.0,gcm.slp,
                     interval.1=c(1950,1999),interval.2=c(2000,2099))

for (im in 1:12) {
  print("catFields")
  XX.t2m <- catFields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- EOF(XX.t2m,plot=FALSE)
  XX.slp <- catFields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- EOF(XX.slp,plot=FALSE)
  print("mixFields:")
  XX <- mixFields(XX.t2m,XX.slp)
  print("mixFields EOF")
  EOF.t2m.slp <- EOF(XX,plot=FALSE)
}


if (do.all) {  # REB 30.09.03

  
##################### ECHAM4/OPYC3: GSDIO ########################

print("### ECHAM4/OPYC3 GSDIO: ###")
gcm.t2m.0 <- retrieve.nc("/home/kareb/data/mpi/mpi-gsdio_t2m.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1860,1990))
gcm.slp.0 <- retrieve.nc("/home/kareb/data/mpi/mpi-gsdio_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1860,1990))

##################### ECHAM4/OPYC3: B2 ########################

print("### ECHAM4/OPYC3 B2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/EH4OPYC_B2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/EH4OPYC_B2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
gcm.t2m <- catFields(gcm.t2m.0,gcm.t2m)
gcm.slp <- catFields(gcm.slp.0,gcm.slp)

for (im in 1:12) {
  XX.t2m <- catFields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- EOF(XX.t2m,plot=FALSE)
  XX.slp <- catFields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- EOF(XX.slp,plot=FALSE)
  XX <- mixFields(XX.t2m,XX.slp)
  EOF.t2m.slp <- EOF(XX,plot=FALSE)
}


##################### ECHAM4/OPYC3: A2 ########################

print("### ECHAM4/OPYC3 A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/EH4OPYC_A2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/EH4OPYC_A2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
for (im in 1:12) {
  XX.t2m <- catFields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- EOF(XX.t2m,plot=FALSE)
  XX.slp <- catFields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- EOF(XX.slp,plot=FALSE)
  XX <- mixFields(XX.t2m,XX.slp)
  EOF.t2m.slp <- EOF(XX,plot=FALSE)
}




##################### ECHAM4/OPYC3: GSDIO ########################

print("### ECHAM4/OPYC3 GSDIO: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/mpi/mpi-gsdio_t2m.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1860,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/mpi/mpi-gsdio_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1860,2099))
for (im in 1:12) {
  XX.t2m <- catFields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- EOF(XX.t2m,plot=FALSE)
  XX.slp <- catFields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- EOF(XX.slp,plot=FALSE)
  XX <- mixFields(XX.t2m,XX.slp)
  EOF.t2m.slp <- EOF(XX,plot=FALSE)
}


  
##################### CSIRO: B2 ########################


print("### CSIRO B2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/CSIRO_B2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1961,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/CSIRO_B2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1961,2099))
for (im in 1:12) {
  XX.t2m <- catFields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- EOF(XX.t2m,plot=FALSE)
  XX.slp <- catFields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- EOF(XX.slp,plot=FALSE)
  XX <- mixFields(XX.t2m,XX.slp)
  EOF.t2m.slp <- EOF(XX,plot=FALSE)
}

##################### CSIRO: A2 ########################

print("### CSIRO A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/CSIRO_A2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1961,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/CSIRO_A2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1961,2099))
for (im in 1:12) {
  XX.t2m <- catFields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- EOF(XX.t2m,plot=FALSE)
  XX.slp <- catFields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- EOF(XX.slp,plot=FALSE)
  XX <- mixFields(XX.t2m,XX.slp)
  EOF.t2m.slp <- EOF(XX,plot=FALSE)
}


##################### CCC: B2 ########################

print("### CCC B2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/CCCma_B2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/CCCma_B2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
for (im in 1:12) {
  XX.t2m <- catFields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- EOF(XX.t2m,plot=FALSE)
  XX.slp <- catFields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- EOF(XX.slp,plot=FALSE)
  XX <- mixFields(XX.t2m,XX.slp)
  EOF.t2m.slp <- EOF(XX,plot=FALSE)
}

##################### CCC: A2 ########################

print("### CCC A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/CCCma_A2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/CCCma_A2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
for (im in 1:12) {
  XX.t2m <- catFields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- EOF(XX.t2m,plot=FALSE)
  XX.slp <- catFields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- EOF(XX.slp,plot=FALSE)
  XX <- mixFields(XX.t2m,XX.slp)
  EOF.t2m.slp <- EOF(XX,plot=FALSE)
}


##################### NCARCSM: A2 ########################

print("### NCARCSM: A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/NCARCSM_A2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/NCARCSM_A2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
for (im in 1:12) {
  XX.t2m <- catFields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- EOF(XX.t2m,plot=FALSE)
  XX.slp <- catFields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- EOF(XX.slp,plot=FALSE)
  XX <- mixFields(XX.t2m,XX.slp)
  EOF.t2m.slp <- EOF(XX,plot=FALSE)
}


##################### NCARPCM: B2 ########################

print("### NCARPCM: B2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/NCARPCM_B2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/NCARPCM_B2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(2001,2099))
for (im in 1:12) {
  XX.t2m <- catFields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- EOF(XX.t2m,plot=FALSE)
  XX.slp <- catFields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- EOF(XX.slp,plot=FALSE)
  XX <- mixFields(XX.t2m,XX.slp)
  EOF.t2m.slp <- EOF(XX,plot=FALSE)
}

##################### NCARPCM: A2 ########################


print("### NCARPCM: A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/NCARPCM_A2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/NCARPCM_A2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
for (im in 1:12) {
  XX.t2m <- catFields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- EOF(XX.t2m,plot=FALSE)
  XX.slp <- catFields(obs.slp,gcm.slp,mon=im)
  EOF.slp <- EOF(XX.slp,plot=FALSE)
  XX <- mixFields(XX.t2m,XX.slp)
  EOF.t2m.slp <- EOF(XX,plot=FALSE)
}

##################### GFDL: B2 ########################

print("### GFDL: A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/GFDL_B2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/GFDL_B2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
for (im in 1:12) {
  XX.t2m <- catFields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- EOF(XX.t2m,plot=FALSE)
#  XX.slp <- catFields(obs.slp,gcm.slp,mon=im)
#  EOF.slp <- EOF(XX.slp,plot=FALSE)
#  XX <- mixFields(XX.t2m,XX.slp)
#  EOF.t2m.slp <- EOF(XX,plot=FALSE)
}

##################### GFDL: A2 ########################

print("### GFDL: A2: ###")
gcm.t2m <- retrieve.nc("/home/kareb/data/ipcc_sres/GFDL_A2_temp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
gcm.slp <- retrieve.nc("/home/kareb/data/ipcc_sres/GFDL_A2_slp.nc",
                       x.rng=x.rng,y.rng=y.rng,t.rng=c(1990,2099))
for (im in 1:12) {
  XX.t2m <- catFields(obs.t2m,gcm.t2m,mon=im)
  EOF.t2m <- EOF(XX.t2m,plot=FALSE)
#  XX.slp <- catFields(obs.slp,gcm.slp,mon=im)
#  EOF.slp <- EOF(XX.slp,plot=FALSE)
#  XX <- mixFields(XX.t2m,XX.slp)
#  EOF.t2m.slp <- EOF(XX,plot=FALSE)
}

} # end do REB 30.09.03
do <- TRUE



}
}

#source("cmp_eofs_gsdio.R")
