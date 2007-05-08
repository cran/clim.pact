# Read the CET from the Internet:
print("Read the CET from the Internet and plot the station data")
a <- readline("[press to proceed]")

cet <- getHadObs()
plotStation(cet)
data(DNMI.slp)

print("plot correlation map")
a <- readline("[press to proceed]")
while (dev.cur()>1) dev.off()

corField(DNMI.slp,cet,mon=1)

# Read on-line station data from the NARP project
print(" Read on-line station data from the NARP project from Internet")
a <- readline("[press to proceed]")
while (dev.cur()>1) dev.off()

NARP <- getnarp()
print(NARP$name)
nuuk <- getnarp("Nuuk")
plotStation(nuuk)

# Construct a 'station object':
print("Construct a 'station object'")
while (dev.cur()>1) dev.off()
a <- readline("[press to proceed]")

rnd <- station.obj(x=matrix(rnorm(100*12),100,12),yy=1901:2000,
                   obs.name="Random",unit="dimensionless",ele=101, mm=NULL,
                   station=NULL,lat=60,lon=0,alt=NULL,
                   location="unspecified",wmo.no=NULL,
                   start=NULL,yy0=NULL,country=NULL,ref=NULL)
# plot (synthetic) random station data
print("plot (synthetic) random station data and correlate with SLP")
while (dev.cur()>1) dev.off()
a <- readline("[press to proceed]")

plotStation(rnd)
corField(DNMI.slp,rnd,mon=1)

# Retrieve SSTs for the North Atlantic.
print("Retrieve SSTs for the North Atlantic.")
a <- readline("[press to proceed]")
while (dev.cur()>1) dev.off()

data(DNMI.sst)
corField(DNMI.sst,nuuk,mon=c(12,1,2))

# Field object handling
print("Load EOF and re-convert to its original form. Use to compute correlations")
a <- readline("[press to proceed]")
while (dev.cur()>1) dev.off()
data(eof.slp)
plotEOF(eof.slp)
SLP <- EOF2field(eof.slp)
corField(SLP,nuuk,mon=1)

print("Plot map of last observations in field and the mean field")
a <- readline("[press to proceed]")
while (dev.cur()>1) dev.off()
mapField(DNMI.sst)
mT2m <- meanField(DNMI.sst)
map(mT2m)

#----------------------------------------------------------------------------------
# ESD for Monthly data:
#----------------------------------------------------------------------------------

# simple ESD:
print("compute EOFs")
a <- readline("[press to proceed]")
while (dev.cur()>1) dev.off()
data(DNMI.t2m)
eof.t2m.1 <- EOF(DNMI.t2m,mon=1)
plotEOF(eof.t2m.1)

# simple ESD:
print("simple ESD")
a <- readline("[press to proceed]")
while (dev.cur()>1) dev.off()
ds.nuuk.1 <- DS(nuuk,eof.t2m.1)

# simple ESD with random values:
print("simple ESD with random values")
a <- readline("[press to proceed]")
while (dev.cur()>1) dev.off()
ds.rnd <- DS(rnd,eof.t2m.1)


