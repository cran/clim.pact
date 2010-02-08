# Overlays coast lines on contour plots/image plots. A mapping tool.
# Assumes that the x-axis and y-axis are given as degrees lon, lat.
#
# Reference: R.E. Benestad et al. (2002),
#            Empirically downscaled temperature scenarios for Svalbard,
#            submitted to Atm. Sci. Lett.
#
#            R.E. Benestad (2001),
#            A comparison between two empirical downscaling strategies,
#            Int. J. Climatology, 1645-1668, vol. 21, DOI 10.1002/joc.703
#
# R.E. Benestad, met.no, Oslo, Norway 16.04.2002
# rasmus.benestad@met.no
#------------------------------------------------------------------------


addland<-function(col="grey50",lwd=1) {

data("addland1",envir = environment())
lines(lon.cont,lat.cont,type="l",col=col,lwd=lwd)
lon.cont<-lon.cont+360
lines(lon.cont,lat.cont,type="l",col=col,lwd=lwd)

}
